cd(@__DIR__)
using Pkg
Pkg.activate(".")
##
using StochasticPrograms
using GLPK
using Pipe
using Random
using Plots
using JuMP
using DataFrames
using CSV
using Cbc

include("samplers.jl")
include("plot_utils.jl")
##
function plot_results_r(sp, pv, w, d; s=1, stage_1=[:gci, :gco], stage_2=[:gci2, :gco2], debug = false)
    plt_sto = plot()
    plt_invest = plot()
    plt = plot() # create separate plot for storage state of charge
    plot!(plt_invest, pv .* value(sp[1, :u_pv]), label="pv")
    plot!(plt_invest, w .* value(sp[1, :u_wind]), label="wind")
    plot!(plt_invest, d, label="demand")
    stor_flow_i = value.(sp[1, :sto_in]).data
    stor_flow_o = value.(sp[1, :sto_out]).data
    t_xi = scenarios(sp)[s].t_xi
    stor_charge = [-sum(stor_flow_i[1:t])+sum(stor_flow_o[1:t]) for t in value.(sp[1, :sto_in]).axes[1]] .+ 0.5*value(sp[1, :u_storage])
    plot!(plt_sto, value.(sp[1, :sto_in]).axes, stor_charge, label="global storage charge")
    if debug
        #print("Maximum value of storage flow = "*string(maximum(value.(sp[1, :sto_in]).data))*"\n")
    end
    for var in stage_1
        plot!(plt, value.(sp[1, var]).axes, value.(sp[1, var]).data, label=string(var))
        if debug
            print("Maximum value of "*string(var)*" = "*string(maximum(value.(sp[1, var]).data))*"\n")
        end
    end
    #=for var in stage_2
        if debug
            print("Maximum value of "*string(var)*" = "*string(maximum(value.(sp[2, var], s).data))*"\n")
        end
        plot!(plt, value.(sp[2, var], s).axes, value.(sp[2, var], s).data, label=string(var)*string(s), linestyle=:dash, linewidth=2)
    end=#
    for var in [:gci2, :gco2]
        plot!(plt, (t_xi+1):(t_xi+3), value.(sp[2, var], s), label=string(var)*string(s), linestyle=:dash, linewidth=2)
    end
    stor_flow_i2 = value.(sp[2, :sto_in2],s)
    stor_flow_o2 = value.(sp[2, :sto_out2],s)    
    stor_charge2 = [-sum(stor_flow_i2[1:t])+sum(stor_flow_o2[1:t]) for t in 1:(reaction_time+1)] .+ stor_charge[t_xi-1]
    plot!(plt_sto, (t_xi-1):(t_xi+reaction_time), vcat(stor_charge[t_xi-1],stor_charge2), label=string("stochastic storage charge")*string(s), linestyle=:dash, linewidth=2)
    if debug
        #print("Maximum value of sto2 = "*string(maximum(value.(sp[2, :sto2], s).data))*"\n")
    end
    display(plot(plt_invest, plt, plt_sto, layout = (3,1)))
    #print(s)
end

## Try adding F to the scenario

@sampler F_sampler_r = begin
    timesteps::UnitRange{Int64}
    F_range::Tuple{Float64, Float64}
    reaction_time::Int64
    #n::Int64
    F_sampler_r(timesteps::UnitRange{Int64}, F_range::Tuple{Float64, Float64}, reaction_time::Int64) = new(timesteps, F_range, reaction_time)
    
    @sample F_scenario begin
        @parameters timesteps F_range
        return F_scenario(rand(timesteps[1:end-reaction_time-1]), rand([-1, 1]), F_range[1]+rand()*F_range[2])
    end
end
##
# Define the time interval
t_s = 1
t_f = 168
offset = 6531
timesteps = t_s:t_f

F_range = (100., 200.)
# Define weather and demand data
pv = CSV.read("../basic_example.csv", DataFrame)[timesteps.+offset, 3]
wind = CSV.read("../basic_example.csv", DataFrame)[timesteps.+offset, 4]
demand = CSV.read("../basic_example.csv", DataFrame)[timesteps.+offset, 2] ./80

##
# Note on signs: s_xi<0 means that energy is requested, s_xi>0 means an additional consumption
energy_model_r = @stochastic_model begin 
    @stage 1 begin
        @parameters begin
            # Expected lifetime of components, years
            asset_lifetime = 10.
            # Costs in Euro/1kWp
            c_pv = 100.
            c_wind = 1000.
            c_storage = 100.
            c_sto_op = 0.00001
            c_i = .03
            c_o = .01
            # Euro
            inv_budget = 500000000.
            # Shiftable demand, kW
            flexible_demand = 0.
        end
        time_scale = asset_lifetime *365*24/length(timesteps)
        # Component units to be invested in, kWp
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
        @objective(model, Min, (u_pv*c_pv+u_wind*c_wind+u_storage*c_storage)/time_scale + c_i*sum(gci) - c_o*sum(gco)+
        c_sto_op*sum(sto_in)+c_sto_op*sum(sto_out))
    end
    @stage 2 begin
        @parameters begin
            c_flex = .5
            #F = 10.
            reaction_time = 3
            c_sto_op = 0.00001
            c_i = .03
            c_o = .01
        end
        @uncertain t_xi s_xi F_xi from F_scenario #t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
        @known(model, u_pv)
        @known(model, u_wind)
        @known(model, u_storage)
        @known(model, gci)
        @known(model, gco)
        @known(model, sto_in)
        @known(model, sto_out)
        @known(model, fl_dem)
        # Post event components
        @recourse(model, gci2[t in 1:reaction_time] >= 0)
        @recourse(model, gco2[t in 1:reaction_time] >= 0)
        @recourse(model, sto_in2[t in 1:(reaction_time+1)] >= 0)
        @recourse(model, sto_out2[t in 1:(reaction_time+1)] >= 0)

        # Storage
        @constraint(model, [t in 1:(reaction_time+1)], 0.5*u_storage - sum(sto_in[t_s:(t_xi-1)])+sum(sto_out[t_s:(t_xi-1)]) <= -sum(sto_in2[t_s:t])+sum(sto_out2[t_s:t]))
        @constraint(model, [t in 1:(reaction_time+1)], -sum(sto_in2[t_s:t])+sum(sto_out2[t_s:t]) <= u_storage - sum(sto_in[t_s:(t_xi-1)])+sum(sto_out[t_s:(t_xi-1)]))

        @constraint(model, -sum(sto_in2)+sum(sto_out2) == -sum(sto_in[t_xi:(t_xi+reaction_time)])+sum(sto_out[t_xi:(t_xi+reaction_time)]))

        #Flexible demand
        @recourse(model, fl_dem2[t in timesteps]>=0)
        @constraint(model, [t in t_s:(t_xi-1)], fl_dem[t] == fl_dem2[t])
        @constraint(model, sum(fl_dem2) == sum(fl_dem))

        # Event energy balance
        # The storage and other fast acting components use the recourse variables here.
        # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
        @constraint(model, gci[t_xi]-gco[t_xi]+u_pv*pv[t_xi]+u_wind*wind[t_xi]-demand[t_xi]-fl_dem2[t_xi] +sto_in2[1]-sto_out2[1]+F_xi*s_xi==0)
        # Post event energy balance
        @constraint(model, [t in (t_xi+1):(t_xi+reaction_time)],
        gci2[t-t_xi]-gco2[t-t_xi]+u_pv*pv[t]+u_wind*wind[t]-demand[t]-fl_dem2[t]+sto_in2[t-t_xi+1]-sto_out2[t-t_xi+1]==0)
        @objective(model, Min, c_i*(sum(gci2)-sum(gci[t_xi:(t_xi+reaction_time)]))-c_o*(sum(gco2)-sum(gco[t_xi:(t_xi+reaction_time)]))+
            c_sto_op*sum(sto_in)+c_sto_op*sum(sto_out)-c_flex*F_xi)
    end
end
## 
reaction_time = 3
f_samp = F_sampler_r(timesteps, F_range, reaction_time)
@timed sp0 = instantiate(energy_model_r, [F_scenario_r(5, 1, 10.)], c_pv = 700., c_wind = 2000., c_i = 0.3, c_o = 0.05, c_storage = 600., flexible_demand = 0., reaction_time = reaction_time, optimizer = GLPK.Optimizer)

##
@timed optimize!(sp0)

od0 = optimal_decision(sp0)
objective_value(sp0)

plot_results_r(sp0, pv, wind, demand, s = 1, stage_1 = [:gci, :gco, :fl_dem], stage_2 = [:gci2, :gco2, :fl_dem2],)

##
@timed sp = instantiate(energy_model_r, f_samp, 100, c_pv = 700., c_wind = 2000., c_i = 0.3, c_o = 0.05, c_storage = 600., flexible_demand = 0., reaction_time = reaction_time, optimizer = GLPK.Optimizer)

@timed optimize!(sp)

od = optimal_decision(sp)
objective_value(sp)

##
#sp = instantiate(energy_model, [simple_scenario(1,0)], F=100., optimizer = Cbc.Optimizer)

##
for s in 1:5
    plot_results_r(sp, pv, wind, demand, s = s, stage_1 = [:gci, :gco, :fl_dem], stage_2 = [:gci2, :gco2, :fl_dem2],)
end
##
# Main result
println("Termination status: $(termination_status(sp))")
println("Objective value: $(objective_value(sp))")
println("Optimal decision: $(optimal_decision(sp))")

##
# Investment result
u_pv = sp[1, :u_pv]
u_wind = sp[1, :u_wind]
u_storage = sp[1, :u_storage]
# First stage
println("value(u_pv) = $(value(u_pv))")
println("value(u_wind) = $(value(u_wind))")
println("value(u_storage) = $(value(u_storage))")

# Second stage
for s in 1:3
    #println("Objective value in scenario $s: $(objective_value(sp, s))")
    println("Optimal recourse in scenario $s: $(optimal_recourse_decision(sp, s))")
    plot_results(sp, pv, wind, demand, s = s)
end

##

scen = scenarios(sp)
print(scen)
##
function test_decision(p, decision, timesteps)

    infeasible_count = 0
    costs = []


    for sign in [-1,1]
        for t in timesteps
            scenario = F_scenario(t,sign,150.)
            try evaluate_decision(p, decision, scenario);
            catch e;
                infeasible_count += 1;
                append!(costs, Inf)
                #append!(scenario_results, [Inf, scenario.t_xi, scenario.s_xi, scen.F_xi]);
            end
            append!(costs, evaluate_decision(p, decision, scenario))
            #append!(scenario_results, [evaluate_decision(p, decision, scenario), scenario.t_xi, scenario.s_xi, scenario.F_xi]);
        end
    end
    N = length(timesteps)
    print(infeasible_count/N/2.)

    scenario_results = hcat(costs, vcat(timesteps, timesteps), -1*vcat(ones(N), ones(N)))
    return scenario_results
end
##

scen_results = test_decision(sp, od, timesteps)

scen_results0 = test_decision(sp0, od0, timesteps)

#av_cost = sum(scen_results[i, 1] for )

print(scen_bad0)

EVPI(sp)

VRP(sp)

VSS(sp)

set_optimizer(sp, SAA.Optimizer)

evaluate_decision(sp, od, simple_scenario(1,1))


t_xi = scenarios(sp)[3].t_xi
