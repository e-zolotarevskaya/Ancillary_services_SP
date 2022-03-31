cd(@__DIR__)
using Pkg
Pkg.activate(".")
##

include("sp_model.jl")
include("plot_utils.jl")
#include("results_handling.jl")
include("evaluation_utils.jl")
##

using DataFrames
using CSV

offset = 6531
timesteps = 1:168

pv = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 3]
wind = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 4]
demand = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps .+ offset, 2]

##

pars = copy(default_es_pars)

##

using Statistics
average_hourly_demand = mean(demand)

pars[:flexible_demand] = average_hourly_demand * 100. # Over the year, we have 100 demand hours that are shiftable freely...
pars[:recovery_time] = 24
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 800.
pars[:c_sto_op] = 0.00001

heatdemand = copy(demand)./300.
heatdemand0 = zeros(length(demand))
es = define_energy_system(pv, wind, demand, heatdemand; p = pars)

##
n = 100
F_max = average_hourly_demand * 0.1 # Have ancillary services equal to 10% of our typical demand
t_max = length(pv) - es.parameters[2].defaults[:recovery_time]
scens = simple_flex_sampler(n, F_max, t_max)

##

using Cbc

sp = instantiate(es, scens, optimizer = Cbc.Optimizer)

# using GLPK

# sp = instantiate(es, scens, optimizer = GLPK.Optimizer)

##

optimize!(sp)

##

println("Objective value: $(objective_value(sp))")

println("PV: $(round(value.(sp[1, :u_pv]); digits = 2)) -- Wind: $(round(value.(sp[1, :u_wind]); digits = 2)) -- Battery: $(round(value.(sp[1, :u_storage]); digits = 2))")
println("Heat storage: $(round(value.(sp[1, :u_heat_storage]); digits = 2))")

# for s in 1:minimum((n, 10))
#     println("Second stage objective value in scenario $s: $(objective_value(sp, s))")
# end

##

plot_results(sp, pv, wind, demand, hd = heatdemand, s = 8, stage_1 = [:gci, :gco], stage_2 = [:gci2, :gco2])

plot_difference(sp, s=8)

##
scen_results = test_decision(sp, timesteps)
scen_results_F = test_decision_variate_F(sp, timesteps, F_step = 500., F_max = 1000.)

plot_flexibility(scen_results, objective_value(sp))
