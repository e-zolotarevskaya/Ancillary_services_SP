cd(@__DIR__)
using Pkg
Pkg.activate(".")

using StochasticPrograms
using GLPK
using Pipe
using Random
using Plots
using JuMP

## Utility functions
function clamped_sin(x)
    return sin(x-8)>0 ? 3. .*sin(x-8) : 0.
end

function wind_timeseries(Length)
    return 4. .+ 2. .*rand(Length)
end

function demand_timeseries(Length)
    return 10. .+2. .*rand(Length)
end

function Theta(x)
    return x>=0 ? 1 : 0
end
##
function plot_results(od, pv, w, d)
    t = Int((length(od)-3)/2+3)
    plt = plot(od[4:t], label = "gci")
    plot!(plt, -od[t+1:end], label = "gco")
    plot!(plt, pv.*od[1], label = "pv")
    plot!(plt, w.*od[2], label = "wind")
    plot!(plt, d, label = "demand")
    display(plt)
end
##
# Global system parameter
c_i = .3
c_o = .1
timesteps = 1:24
time_scale = 25. *365*24/length(timesteps)
c_pv = 1000.
c_wind = 1000.
c_storage = 400.
c_flex = .5
F = 3
storage_scale = 20.
pv = [clamped_sin(x/length(timesteps) * 2. * pi) for x in timesteps]
wind = wind_timeseries(length(timesteps))
demand = demand_timeseries(length(timesteps))
##
function scenarios(n, timesteps)
    return [@scenario t_xi = rand(1:length(timesteps)) s_xi = rand([0,1]) probability = 1. /n for i in 1:n]
end
##
#= Note on signs: s_xi<0 means that energy is requested, s_xi>0 means an additional consumption
To simplify the problem and avoid using inequalities in constraints, let's first study only s_xi<0
Then buy-back B>0, therefore it's associated with price c_i
=#
@stochastic_model em begin
    @stage 1 begin
        @decision(em, u_pv>=0)
        @decision(em, u_wind>=0)
        @decision(em, u_storage>=0)
        @decision(em, gci[t in timesteps]>=0)
        @decision(em, gco[t in timesteps]>=0)
        @objective(em, Min, u_pv*c_pv/time_scale+u_wind*c_wind/time_scale+u_storage*c_storage/time_scale+c_i*sum(gci)-c_o*sum(gco))
    end
    @stage 2 begin
        @uncertain t_xi s_xi #t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
        @recourse(em, storage[t in timesteps])
        @recourse(em, B[t in timesteps])
        @constraint(em, [t in timesteps], s_xi*B[t]<=0)
        @constraint(em, [t in 1:t_xi], B[t]==0)
        @constraint(em, [t in 1:(t_xi-1)], 
        storage[t] == gco[t]-gci[t]-u_wind*wind[t]-c_pv*pv[t]+demand[t])
        @constraint(em, 
        storage[t_xi]+F*s_xi+B[t_xi] == gco[t_xi]-gci[t_xi]-u_wind*wind[t_xi]-c_pv*pv[t_xi]+demand[t_xi])
        @constraint(em, [t in (t_xi+1):length(timesteps)], 
        storage[t]+B[t] == gco[t]-gci[t]-u_wind*wind[t]-c_pv*pv[t]+demand[t])
        @objective(em, Min, c_i*sum(B)+u_pv*c_pv/time_scale+u_wind*c_wind/time_scale+u_storage*c_storage/time_scale+c_i*sum(gci)-c_o*sum(gco))
        #@constraint(em, b[1] == 0.5*c_b*e_c) # initial charge
        @constraint(em, [t in timesteps], -0.5*u_storage*storage_scale <= sum(storage[1:t]))
        @constraint(em, [t in timesteps], sum(storage[1:t]) <= 0.5*u_storage*storage_scale)
    end
end

## 
xi = scenarios(10, timesteps)
xi_1 = @scenario t_xi = 2 s_xi = -1 probability = 0.5

xi_2 = @scenario t_xi = 10 s_xi = -1 probability = 0.5
sp = instantiate(em, [x for x in xi], optimizer = GLPK.Optimizer)
##
optimize!(sp)
objective_value(sp)

od = optimal_decision(sp)
##
#=xi = scenarios(2, 24)
sp = instantiate(em, [x for x in xi], optimizer = GLPK.Optimizer)=#
plot_results(od, pv, wind, demand)

##
# Main result
println("Termination status: $(termination_status(sp))")
println("Objective value: $(objective_value(sp))")
println("Optimal decision: $(optimal_decision(sp))")

# First stage
println("value(u_pv) = $(value(u_pv))")
println("reduced_cost(x) = $(reduced_cost(x))")

# Scenario 1
# Second stage
println("value(y, 1) = $(value(y, 1))")
println("reduced_cost(y, 1) = $(reduced_cost(y, 1))")
println("dual(con, 1) = $(dual(con, 1))")
println("Objective value in scenario 1: $(objective_value(sp, 1))")
println("Optimal recourse in scenario 1: $(optimal_recourse_decision(sp, 1))")

# Scenario 2
println("value(y, 2) = $(value(y, 2))")
println("reduced_cost(y, 2) = $(reduced_cost(y, 2))")
println("dual(con, 2) = $(dual(con, 2))")
println("Objective value in scenario 2: $(objective_value(sp, 2))")
println("Optimal recourse in scenario 2: $(optimal_recourse_decision(sp, 2))")
## 
#=
model em begin

        @variable(em, c_pv>=0)
        @variable(em, c_w>=0)
        @variable(em, c_b>=0)
        @variable(em, gci[t in timesteps]>=0)
        @variable(em, gco[t in timesteps]>=0)
        @variable(em, b[t in timesteps])
        @variable(em, B[t in timesteps])
        @constraint(em, [t in timesteps], s_xi*B[t]<=0)
        @constraint(em, [t in 1:t_xi], B[t]==0)
        @constraint(em, [t in 1:(t_xi-1)], 
        b[t] == gco[t]-gci[t]-c_w*w[t]-c_pv*pv[t]+d[t])
        @constraint(em, 
        b[t_xi]+F*s_xi+B[t_xi] == gco[t_xi]-gci[t_xi]-c_w*w[t_xi]-c_pv*pv[t_xi]+d[t_xi])
        @constraint(em, [t in (t_xi+1):length(timesteps)], 
        b[t]+B[t] == gco[t]-gci[t]-c_w*w[t]-c_pv*pv[t]+d[t])
        @objective(em, Min, c_i*sum(B)+c_pv/n_y+c_w/n_y+c_b/n_y+c_i*sum(gci)-c_o*sum(gco))
        #@constraint(em, b[1] == 0.5*c_b*e_c) # initial charge
        @constraint(em, [t in timesteps], -0.5*c_b*e_c <= sum(b[1:t]))
        @constraint(em, [t in timesteps], sum(b[1:t]) <= 0.5*c_b*e_c)
end=#