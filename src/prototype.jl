cd(@__DIR__)
using Pkg
Pkg.activate(".")

using StochasticPrograms
using GLPK
using Pipe
using Random
using Plots
using JuMP
using DataFrames
using CSV

function plot_res2(sp, pv, w, d; scenarios=[1], stage_1=[:gci, :gco, :storage], stage_2=[:gci2, :gco2])
    plt = plot(value(sp[1, :u_pv]) .* pv, label="pv")
    plot!(plt, value(sp[1, :u_wind]) .*w, label="wind")
    plot!(plt, d, label="demand")
    for var in stage_1
        plot!(plt, value.(sp[1, var]).axes, value.(sp[1, var]).data, label=string(var))
    end
    for s in scenarios
        for var in stage_2
            plot!(plt, value.(sp[2, var], s).axes, value.(sp[2, var], s).data, label=string(var))
        end
    end
    plt
end

##
function plot_results(sp, pv, w, d; scen = nothing)
    od = optimal_decision(sp)
    t = Int((length(od)-3)/2+3)
    gci = od[4:t]
    gco = od[t+1:end]
    gc = [gci[i]==0 ? -gco[i] : gci[i] for i in 1:length(gci)]
    plt = plot(gc, label = "grid connection")
    if od[1] != 0
        plot!(plt, pv.*od[1], label = "pv")
    end
    if od[2] != 0 
        plot!(plt, w.*od[2], label = "wind")
    end
    plot!(plt, -d, label = "demand", xlabel = "t (h)", ylabel = "P (kW)")
    if !isnothing(scen)
    println("Objective value in scenario $scen: $(objective_value(sp, scen))")
    #t = Int((length(od)-3)/2+3)
        ord = optimal_recourse_decision(sp, scen)
        stor_flow = ord[1:t]
        stor = [sum(stor_flow[1:i]) for i in 1:t]
        stor.+=0.5*od[3]
        buyback = ord[t:end]
        plot!(plt, buyback, label = "buy-back")
        plt_st = plot(stor, label = "storage charge")
        plt_general = plot(plt, plt_st, layout = (2, 1))
        display(plt_general)
    else display(plt)
    end
end
##
# Global system parameters
c_i = .03
c_o = .01
t_f = 48
t_s = 1
timesteps = t_s:t_f
time_scale = 10. *365*24/length(timesteps)
c_pv = 1000.
c_wind = 1000.
c_storage = 400.
c_flex = .5
F = 100.
storage_scale = 1.
##
# Define weather and demand data
pv = CSV.read("../data/pv_Halle18.csv", DataFrame)[timesteps, 1]
wind = CSV.read("../data/wind_Karholz.csv", DataFrame)[timesteps, 1]
demand = CSV.read("../data/demand_Industriepark.csv", DataFrame)[timesteps, 1]
#=pv = CSV.read("../data/basic_example.csv", DataFrame)[timesteps, 3]
wind = CSV.read("../data/basic_example.csv", DataFrame)[timesteps, 4]
demand = CSV.read("../data/basic_example.csv", DataFrame)[timesteps, 2]
=#
##
function scenarios(n, timesteps)
    return [@scenario t_xi = rand(1:length(timesteps)-2) s_xi = rand([-1]) probability = 1. /n for i in 1:n]
end
##
#= Note on signs: s_xi<0 means that energy is requested, s_xi>0 means an additional consumption
To simplify the problem and avoid using inequalities in constraints, let's first study only s_xi<0
Then buy-back B>0, therefore it's associated with price c_i
=#
@stochastic_model em begin
    @stage 1 begin
        # Investments:
        @decision(em, u_pv>=0)
        @decision(em, u_wind>=0)
        @decision(em, u_storage>=0)
        # Grid connection
        @decision(em, gci[t in timesteps]>=0)
        @decision(em, gco[t in timesteps]>=0)
        # Storage model
        @decision(em, storage[t in timesteps])
        @constraint(em, [t in timesteps], -0.5*u_storage*storage_scale <= sum(storage[1:t]))
        @constraint(em, [t in timesteps], sum(storage[1:t]) <= 0.5*u_storage*storage_scale)
        @constraint(em, sum(storage) == 0)
        # Energy balance
        @constraint(em, [t in timesteps], 
        gci[t]-gco[t]+u_pv*pv[t]+u_wind*wind[t]-demand[t]+storage[t]==0) # Energy balance
        # Investment costs
        @objective(em, Min, u_pv*c_pv/time_scale+u_wind*c_wind/time_scale+u_storage*c_storage/time_scale + c_i*sum(gci) - c_o*sum(gco))
    end
    @stage 2 begin
        @uncertain t_xi s_xi #t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
        @known(em, u_pv)
        @known(em, u_wind)
        @known(em, u_storage)
        @known(em, gci)
        @known(em, gco)
        @known(em, storage)
        # Post event components
        @recourse(em, gci2[t in t_xi+1:t_f])
        @recourse(em, gco2[t in t_xi+1:t_f])
        @recourse(em, sto2[t in timesteps])
        @constraint(em, [t in timesteps], -0.5*u_storage*storage_scale <= sum(sto2[1:t]))
        @constraint(em, [t in timesteps], sum(sto2[1:t]) <= 0.5*u_storage*storage_scale)
        @constraint(em, sum(sto2) == 0)

        # Put storage at time of event into same state
        @constraint(em, sum(storage[t_s:t_xi-1]) == sum(sto2[t_s:t_xi-1]))
        
        #@constraint(em, [t in timesteps], 0 <= sum(sto2[t_xi:t_f] +0.5*u_storage*storage_scale + sum(storage[1:t_xi])))
        #@constraint(em, [t in timesteps], sum(sto2[t_xi:t_f]) <= 0.5*u_storage*storage_scale)
        #@constraint(em, sum(sto2) == sum(storage[1:t_xi]))
        
        # Event energy balance
        # The storage and other fast acting components use the recourse variables here.
        # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
        @constraint(em, gci[t_xi]-gco[t_xi]+u_pv*pv[t_xi]+u_wind*wind[t_xi]-demand[t_xi]+sto2[t_xi]+F*s_xi==0)
        # Post event energy balance
        @constraint(em, [t in (t_xi+1):t_f],
        gci2[t]-gco2[t]+u_pv*pv[t]+u_wind*wind[t]-demand[t]+sto2[t]==0)
        @objective(em, Min, c_i*sum(gci2[t_xi+1:t_f])-c_o*sum(gco2[t_xi+1:t_f]) - (c_i*sum(gci[t_xi+1:t_f])-c_o*sum(gco[t_xi+1:t_f])) )
    end
end
## 
#xi = scenarios(10, timesteps)
xi_1 = @scenario t_xi = 3 s_xi = -1 probability = 0.5
xi_2 = @scenario t_xi = 3 s_xi = -1 probability = 0.5
xi = scenarios(10, timesteps)
sp = instantiate(em, xi, optimizer = GLPK.Optimizer)
##
optimize!(sp)

od = optimal_decision(sp)
objective_value(sp)
##
plot_results(sp, pv, wind, demand)

##
# Main result
println("Termination status: $(termination_status(sp))")
println("Objective value: $(objective_value(sp))")
println("Optimal decision: $(optimal_decision(sp))")

##
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
    plot_results(sp, pv, wind, demand, scen = s)
end

##

plot_results(sp, pv, wind, demand, scen = 1)

ord = optimal_recourse_decision(sp, 1)

DeterministicEquivalent(sp)