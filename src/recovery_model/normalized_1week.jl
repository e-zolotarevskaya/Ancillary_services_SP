cd(@__DIR__)
using Pkg
Pkg.activate(".")
##

include("sp_model.jl")
include("plot_utils.jl")

##

using DataFrames
using CSV

timesteps = 1:168

pv = CSV.read("timeseries/basic_example_normalized.csv", DataFrame)[timesteps, 3]
wind = CSV.read("timeseries/basic_example_normalized.csv", DataFrame)[timesteps, 4]
demand = CSV.read("timeseries/basic_example_normalized.csv", DataFrame)[timesteps, 2]

##

pars = copy(default_es_pars)

##

pars[:flexible_demand] = 1000.
es = define_energy_system(pv, wind, demand; p = pars)

##
n = 40
F_max = 100.
t_max = length(pv) - es.parameters[2].defaults[:recovery_time]
scen = simple_flex_sampler(n, F_max, t_max)

##

using GLPK

sp = instantiate(es, scen, optimizer = GLPK.Optimizer)

##

optimize!(sp)

# od = optimal_decision(sp)
println("Objective value: $(objective_value(sp))")

for s in 1:n
    println("Second stage objective value in scenario $s: $(objective_value(sp, s))")
end

##

plot_results(sp, pv, wind, demand, s = 8)


## Run the optimization with the LShaped optimizer

spL = instantiate(es, scen, optimizer = LShaped.Optimizer)

set_optimizer_attribute(spL, MasterOptimizer(), GLPK.Optimizer)
set_optimizer_attribute(spL, SubProblemOptimizer(), GLPK.Optimizer)
set_optimizer_attribute(spL, FeasibilityStrategy(), FeasibilityCuts())

optimize!(spL)

println("Objective value: $(objective_value(spL))")

for s in 1:n
    println("Second stage objective value in scenario $s: $(objective_value(spL, s))")
end

plot_results(spL, pv, wind, demand, s = 8)
