cd(@__DIR__)
using Pkg
Pkg.activate(".")
##
using StochasticPrograms
using GLPK
using Pipe
using Random
using JuMP
using DataFrames
using CSV
using TimerOutputs
##
@define_scenario simple_scenario = begin
    t_xi::Int64
    s_xi::Int64
    @zero begin
        return simple_scenario(0, 0)
    end
    @expectation begin
        t_xi = Int(floor(sum([probability(s)*s.t_xi for s in scenarios])))
        s_xi = Int(floor(sum([probability(s)*s.s_xi for s in scenarios])))
        return simple_scenario(t_xi, s_xi)
    end
end

@sampler simple_sampler = begin
    timesteps::UnitRange{Int64}
    #n::Int64
    simple_sampler(timesteps::UnitRange{Int64}) = new(timesteps)
    
    @sample simple_scenario begin
        @parameters timesteps
        return simple_scenario(rand(timesteps[1:end-2]), rand([-1, 1]))
    end
end
##


##
# Define weather and demand data
function create_sampler(timesteps)
    pv = CSV.read("../basic_example.csv", DataFrame)[timesteps, 3]
    wind = CSV.read("../basic_example.csv", DataFrame)[timesteps, 4]
    demand = CSV.read("../basic_example.csv", DataFrame)[timesteps, 2]
    s = simple_sampler(timesteps)
    return pv, wind, demand, s
end

##

# Note on signs: s_xi<0 means that energy is requested, s_xi>0 means an additional consumption
energy_model = @stochastic_model begin 
    @stage 1 begin
        @parameters begin
            asset_lifetime = 10.
            c_pv = 100.
            c_wind = 1000.
            c_storage = 100.
            inv_budget = 500000000.
            flexible_demand = 0.
        end
        time_scale = asset_lifetime *365*24/length(timesteps)
        # Investments:
        @decision(model, u_pv >= 0)
        @decision(model, u_wind >= 0)
        @decision(model, u_storage >= 0)
        @constraint(model, c_pv*u_pv+c_wind*u_wind+c_storage*u_storage<=inv_budget)
        # Grid connection
        @decision(model, gci[t in timesteps] >= 0)
        @decision(model, gco[t in timesteps] >= 0)
        # Flexible demand
        @decision(model, fl_dem[t in timesteps] >= 0)
        @constraint(model, sum(fl_dem) == flexible_demand)
        # Storage model
        @decision(model, sto_in[t in timesteps] >= 0) #into the bus from storage
        @decision(model, sto_out[t in timesteps] >= 0)
        @constraint(model, [t in timesteps], -0.5*u_storage <= -sum(sto_in[t_s:t])+sum(sto_out[t_s:t]))
        @constraint(model, [t in timesteps], -sum(sto_in[t_s:t])+sum(sto_out[t_s:t]) <= 0.5*u_storage)
        @constraint(model, -sum(sto_in)+sum(sto_out) == 0)
        # Energy balance
        @constraint(model, [t in timesteps], 
        gci[t]-gco[t]+u_pv*pv[t]+u_wind*wind[t]-demand[t]-fl_dem[t]+sto_in[t]-sto_out[t]==0)
        # Investment costs
        @objective(model, Min, (u_pv*c_pv+u_wind*c_wind+u_storage*c_storage)/time_scale)
    end
    @stage 2 begin
        @parameters begin
            c_i = .03
            c_o = .01
            c_flex = .5
            #F = 10.
            c_sto_op = 0.00001
            F_xi = 100.
        end
        @uncertain t_xi s_xi from simple_scenario #t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
        @known(model, u_pv)
        @known(model, u_wind)
        @known(model, u_storage)
        @known(model, gci)
        @known(model, gco)
        @known(model, sto_in)
        @known(model, sto_out)
        @known(model, fl_dem)
        # Post event components
        @recourse(model, gci2[t in t_s:t_f] >= 0)
        @recourse(model, gco2[t in t_s:t_f] >= 0)
        @recourse(model, sto_in2[t in timesteps] >= 0)
        @recourse(model, sto_out2[t in timesteps] >= 0)

        @constraint(model, [t in t_s:t_xi], gci[t] == gci2[t])
        @constraint(model, [t in t_s:t_xi], gco[t] == gco2[t])

        @constraint(model, [t in timesteps], -0.5*u_storage <= -sum(sto_in2[t_s:t])+sum(sto_out2[t_s:t]))
        @constraint(model, [t in timesteps], -sum(sto_in2[t_s:t])+sum(sto_out2[t_s:t]) <= 0.5*u_storage)

        @constraint(model, -sum(sto_in2)+sum(sto_out2) == 0)

        #Flexible demand
        @recourse(model, fl_dem2[t in timesteps]>=0)
        @constraint(model, [t in t_s:(t_xi-1)], fl_dem[t] == fl_dem2[t])
        @constraint(model, sum(fl_dem2) == sum(fl_dem))
        # Put storage at time of event into same state
        @constraint(model, [t in t_s:(t_xi-1)], sto_in[t] == sto_in2[t])
        @constraint(model, [t in t_s:(t_xi-1)], sto_out[t] == sto_out2[t])

        # Event energy balance
        # The storage and other fast acting components use the recourse variables here.
        # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
        @constraint(model, gci[t_xi]-gco[t_xi]+u_pv*pv[t_xi]+u_wind*wind[t_xi]-demand[t_xi]-fl_dem2[t_xi] +sto_in2[t_xi]-sto_out2[t_xi]+F_xi*s_xi==0)
        # Post event energy balance
        @constraint(model, [t in (t_xi+1):t_f],
        gci2[t]-gco2[t]+u_pv*pv[t]+u_wind*wind[t]-demand[t]-fl_dem2[t]+sto_in2[t]-sto_out2[t]==0)
        @objective(model, Min, c_i*sum(gci2)-c_o*sum(gco2)+c_sto_op*sum(sto_in2) + c_sto_op*sum(sto_out2) - c_flex*F_xi)
    end
end
## 

# First try short time intervals with different number of scenarios:
t_s = 1
t_f = 48
timesteps = t_s:t_f
pv, wind, demand, s = create_sampler(timesteps)

@timed sp0 = instantiate(energy_model, s, 1, c_pv = 200., flexible_demand = 0., optimizer = GLPK.Optimizer)
@timed sp = instantiate(energy_model, s, 100, c_pv = 200., flexible_demand = 0., optimizer = GLPK.Optimizer)

# 1 scenario - 0.03 s, 100 scenarios - 0.75 s

# Try with a longer time interval:
t_s = 1
t_f = 480
timesteps = t_s:t_f
pv, wind, demand, s = create_sampler(timesteps)

@timed sp0 = instantiate(energy_model, s, 1, c_pv = 200., flexible_demand = 0., optimizer = GLPK.Optimizer)
@timed sp = instantiate(energy_model, s, 100, c_pv = 200., flexible_demand = 0., optimizer = GLPK.Optimizer)

# 1 scenario - 0.8 s, 100 scenarios - 21.5 s

# Even longer:
t_s = 1
t_f = 4800
timesteps = t_s:t_f
pv, wind, demand, s = create_sampler(timesteps)

@timed sp0 = instantiate(energy_model, s, 1, c_pv = 200., flexible_demand = 0., optimizer = GLPK.Optimizer)
@timed sp = instantiate(energy_model, s, 100, c_pv = 200., flexible_demand = 0., optimizer = GLPK.Optimizer)

# 1 scenario - fails.


##
optimize!(sp0)

od0 = optimal_decision(sp0)
objective_value(sp0)

plot_results(sp0, pv, wind, demand, s = 1, stage_1 = [:gci, :gco, :fl_dem], stage_2 = [:gci2, :gco2, :fl_dem2],)


##
optimize!(sp)

od = optimal_decision(sp)
objective_value(sp)

##
#sp = instantiate(energy_model, [simple_scenario(1,0)], F=100., optimizer = Cbc.Optimizer)
