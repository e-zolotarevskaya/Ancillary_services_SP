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

function plot_results(sp, pv, w, d; scenarios=[1], stage_1=[:gci, :gco], stage_2=[:gci2, :gco2], debug = false)
    plt_sto = plot() # create separate plot for storage state of charge
    plt = plot(pv .* value(sp[1, :u_pv]), label="pv")
    plot!(plt, w .* value(sp[1, :u_wind]), label="wind")
    plot!(plt, -d, label="demand")
    plot!(plt_sto, value.(sp[1, :storage]).axes, value.(sp[1, :storage]).data, label="global storage charge")
    if debug
        print("Maximum value of storage flow = "*string(maximum(value.(sp[1, :storage]).data))*"\n")
    end
    for var in stage_1
        plot!(plt, value.(sp[1, var]).axes, value.(sp[1, var]).data, label=string(var))
        if debug
            print("Maximum value of "*string(var)*" = "*string(maximum(value.(sp[1, var]).data))*"\n")
        end
    end
    for s in scenarios
        for var in stage_2
            plot!(plt, value.(sp[2, var], s).axes, value.(sp[2, var], s).data, label=string(var)*string(s))
            if debug
                print("Maximum value of "*string(var)*" = "*string(maximum(value.(sp[2, var], s).data))*"\n")
            end
        end
        plot!(plt_sto, value.(sp[2, :sto2], s).axes, value.(sp[2, :sto2], s).data, label=string("stochastic storage charge")*string(s))
        if debug
            print("Maximum value of sto2 = "*string(maximum(value.(sp[2, :sto2], s).data))*"\n")
        end
    end
    display(plot(plt, plt_sto, layout = (2,1)))
end


##
# Global system parameters
c_i = .03
c_o = .01
t_s = 1
t_f = 48
offset = 2400
timesteps = t_s:t_f
time_scale = 10. *365*24/length(timesteps)
c_pv = 100.
c_wind = 1000.
c_storage = 100.
c_flex = .5
F = 10.
storage_scale = 1.
##
# Define weather and demand data
#=pv = CSV.read("../data/pv_Halle18.csv", DataFrame)[timesteps, 1]
wind = CSV.read("../data/wind_Karholz.csv", DataFrame)[timesteps, 1]
demand = CSV.read("../data/demand_Industriepark.csv", DataFrame)[timesteps, 1]=#
pv = CSV.read("../data/basic_example.csv", DataFrame)[timesteps.+offset, 3]
wind = CSV.read("../data/basic_example.csv", DataFrame)[timesteps.+offset, 4]
demand = CSV.read("../data/basic_example.csv", DataFrame)[timesteps.+offset, 2]

##
function scenarios(n, timesteps)
    return [@scenario t_xi = rand(timesteps[1:end-2]) s_xi = rand([-1, 1]) probability = 1. /n for i in 1:n]
end
##
# Note on signs: s_xi<0 means that energy is requested, s_xi>0 means an additional consumption
@stochastic_model em begin
    @stage 1 begin
        # Investments:
        @decision(em, 100000 >= u_pv >= 0)
        @decision(em, 100000 >= u_wind >= 0)
        @decision(em, 100000 >= u_storage >= 0)
        # Grid connection
        @decision(em, gci[t in timesteps] >= 0)
        @decision(em, gco[t in timesteps] >= 0)
        # Storage model
        @decision(em, storage[t in timesteps])
        @constraint(em, [t in timesteps], -0.5*u_storage*storage_scale <= sum(storage[t_s:t]))
        @constraint(em, [t in timesteps], sum(storage[t_s:t]) <= 0.5*u_storage*storage_scale)
        @constraint(em, sum(storage) == 0)
        # Energy balance
        @constraint(em, [t in timesteps], 
        gci[t]-gco[t]+u_pv*pv[t]+u_wind*wind[t]-demand[t]+storage[t]==0) # Energy balance
        # Investment costs
        @objective(em, Min, u_pv*c_pv/time_scale+u_wind*c_wind/time_scale+u_storage*c_storage/time_scale+c_i*sum(gci)-c_o*sum(gco))
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
        @recourse(em, gci2[t in t_xi+1:t_f] >= 0)
        @recourse(em, gco2[t in t_xi+1:t_f] >= 0)
        @recourse(em, sto2[t in timesteps])
        @constraint(em, [t in timesteps], -0.5*u_storage*storage_scale <= sum(sto2[t_s:t]))
        @constraint(em, [t in timesteps], sum(sto2[t_s:t]) <= 0.5*u_storage*storage_scale)
        @constraint(em, sum(sto2) == 0)

        # Put storage at time of event into same state
        @constraint(em, sum(storage[t_s:t_xi-1]) == sum(sto2[t_s:t_xi-1]))
        
        #=
        @recourse(em, sto2[t in t_xi+1:t_f])
        @constraint(em, [t in 1:t_xi-1], sto2[t] == storage[t])
        @constraint(em, [t in timesteps], 0 <= sum(sto2[t_xi:t_f]) + sum(storage[1:t_xi]))
        @constraint(em, [t in timesteps], sum(sto2[t_xi:t_f]) + sum(storage[1:t_xi]) <= u_storage*storage_scale)
        @constraint(em, sum(sto2) == - sum(storage[1:t_xi]))
        =#
        # Event energy balance
        # The storage and other fast acting components use the recourse variables here.
        # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
        @constraint(em, gci[t_xi]-gco[t_xi]+u_pv*pv[t_xi]+u_wind*wind[t_xi]-demand[t_xi]+sto2[t_xi]+F*s_xi==0)
        # Post event energy balance
        @constraint(em, [t in (t_xi+1):t_f],
        gci2[t]-gco2[t]+u_pv*pv[t]+u_wind*wind[t]-demand[t]+sto2[t]==0)
        @objective(em, Min, c_i*sum(gci2[t_xi+1:t_f])-c_o*sum(gco2[t_xi+1:t_f]) - (c_i*sum(gci[t_xi+1:t_f])-c_o*sum(gco[t_xi+1:t_f])))
    end
end
## 
#xi_1 = @scenario t_xi = 5 s_xi = -1 probability = 0.5
#xi_2 = @scenario t_xi = 3 s_xi = -1 probability = 0.5 #Int(floor((t_f+t_s)/2))
xi = scenarios(10, t_s:t_f)
sp = instantiate(em, xi, optimizer = GLPK.Optimizer)
##
optimize!(sp)

od = optimal_decision(sp)
objective_value(sp)
##
plot_results(sp, pv, wind, demand, debug = true)

for t in timesteps
    xi_t = @scenario t_xi = t s_xi = -1 probability = 1.
    sp_t = instantiate(em, [xi_t], optimizer = GLPK.Optimizer)
    optimize!(sp_t)
    od = optimal_decision(sp_t)
    print(objective_value(sp_t))
    print(", t_xi = $t \n")
end
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
    plot_results(sp, pv, wind, demand, scenarios = [s])
end

##

plot_results(sp, pv, wind, demand, scenarios = 1:2)