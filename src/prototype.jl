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
##

include("samplers.jl")
include("plot_utils.jl")

##
# Define the time interval
t_s = 1
t_f = 168
offset = 6531
timesteps = t_s:t_f

F_range = (100., 200.)
f_samp = F_sampler(timesteps, F_range)
# Define weather and demand data
pv = CSV.read("../basic_example.csv", DataFrame)[timesteps.+offset, 3]
wind = CSV.read("../basic_example.csv", DataFrame)[timesteps.+offset, 4]
demand = CSV.read("../basic_example.csv", DataFrame)[timesteps.+offset, 2] ./80
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
#c_i = 0.5, c_o = 0.04, c_wind = 10., c_pv = 10., c_storage = 1. # very expensive grid connection, cheap storage
#xi = [simple_scenario(5,1), simple_scenario(5,-1)]

sp0 = instantiate(energy_model, [F_scenario(t_f, 0, 0.)], c_pv = 700., c_wind = 2000., c_i = 0.3, c_o = 0.05, c_storage = 600., flexible_demand = 0., optimizer = GLPK.Optimizer)

optimize!(sp0)

od0 = optimal_decision(sp0)
ov0 = objective_value(sp0)


plot_results(sp0, pv, wind, demand, s = 1, stage_1 = [:gci, :gco, :fl_dem], stage_2 = [:gci2, :gco2, :fl_dem2],)

##
sp = instantiate(energy_model, f_samp, 100,  c_pv = 700., c_wind = 2000., c_i = 0.3, c_o = 0.05, c_storage = 600., flexible_demand = 0., optimizer = GLPK.Optimizer)
optimize!(sp)

od = optimal_decision(sp)
ov = objective_value(sp)

##
#sp = instantiate(energy_model, [simple_scenario(1,0)], F=100., optimizer = Cbc.Optimizer)

##
for s in 1:5
    plot_results(sp, pv, wind, demand, s = s, stage_1 = [:gci, :gco, :fl_dem], stage_2 = [:gci2, :gco2, :fl_dem2],)
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
    println("Objective value in scenario $s: $(objective_value(sp, s))")
    #println("Optimal recourse in scenario $s: $(optimal_recourse_decision(sp, s))")
    plot_results(sp, pv, wind, demand, s = s)
end

##

scen = scenarios(sp)
print(scen)
##
function test_decision(p, decision, timesteps)
    infeasible_count = 0
    costs = []
    #scenario_results = []
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

    scenario_results = hcat(costs, vcat(timesteps, timesteps), vcat(-1*ones(N), ones(N)))
    return scenario_results
end

function test_decision_variate_F(p, decision, timesteps; F_step = 50., F_max = 300.)
    costs = zeros(2*length(timesteps))
    F_potential = zeros(2*length(timesteps))
    #scenario_results = []
    for Sign in [-1,1]
        for t in timesteps
            i = t
            if Sign == 1
                i = t + length(timesteps)
            end
            for F in F_step:F_step:F_max
                scenario = F_scenario(t,Sign,F*1.)
                c = evaluate_decision(p, decision, scenario)
                if c!= Inf
                    F_potential[i] = F*Sign
                    costs[i] = c
                else
                    if costs[i] == 0
                        costs[i] = Inf
                        F_potential[i] = F*Sign
                        break
                    end
                end
            end
        end
    end
    print(F_potential)
    scenario_results = hcat(costs, vcat(timesteps, timesteps), F_potential)
    return scenario_results
end

##

scen_results = test_decision_variate_F(sp, od, timesteps)

scen_results0 = test_decision_variate_F(sp0, od0, timesteps)

positive_flexibility, negative_flexibility, positive_potential, negative_potential = plot_flexibility(scen_results0, ov0)

plot_flexibility_cost(scen_results0, ov0)


plot(positive_potential)
plot!(negative_potential)



plot_flexibility(scen_results, ov0)

plot_flexibility(scen_results0, ov0)

##
scen_results = test_decision_variate_F(sp, od, timesteps, F_step = 300., F_max = 3000.)

plt = plot()
plot!(plt, positive_flexibility./positive_potential, label = "price of positive flexibility")
plot!(plt, negative_flexibility./negative_potential, label = "price of negative flexibility", legend = :topleft)
plot!(twinx(), positive_potential, fillrange = 0, fillalpha = 0.35, label = "positive flexibility potential", ylimits = (-300., 300.))
plot!(twinx(), negative_potential, color = :red, fillrange = 0, fillalpha = 0.35, label = "negative flexibility potential", ylimits = (-300., 300.))
#display(plot(plt_cost, plt_pot, layout = (2,1)))
display(plt)