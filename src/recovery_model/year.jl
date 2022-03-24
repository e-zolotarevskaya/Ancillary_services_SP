cd(@__DIR__)
using Pkg
Pkg.activate(".")
##

include("sp_model.jl")
include("plot_utils.jl")

##

using DataFrames
using CSV

timesteps = 1:365*24

pv = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps, 3]
wind = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps, 4]
demand = CSV.read("timeseries/basic_example.csv", DataFrame)[timesteps, 2]

##

pars = copy(default_es_pars)

##

pars[:flexible_demand] = 1000.
pars[:recovery_time] = 24
pars[:c_storage] = 200.
pars[:c_pv] = 300.
pars[:c_wind] = 800.
pars[:c_sto_op] = 0.

es = define_energy_system(pv, wind, demand; p = pars)

##
n = 100
F_max = 100.
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

println("PV: $(value.(sp[1, :u_pv])) -- Wind: $(value.(sp[1, :u_wind])) -- Battery: $(value.(sp[1, :u_storage]))")

# for s in 1:minimum((n, 10))
#     println("Second stage objective value in scenario $s: $(objective_value(sp, s))")
# end

##

plot_results(sp, pv, wind, demand, s = 8)

##