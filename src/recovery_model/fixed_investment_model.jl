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

pars[:recovery_time] = 24
pars[:c_storage] = 100.
pars[:c_pv] = 300.
pars[:c_wind] = 800.
pars[:c_sto_op] = 0.00001

heatdemand = copy(demand)./300.
heatdemand0 = zeros(length(demand))
es = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)


##

using Cbc

sp = instantiate(es, no_flex_pseudo_sampler(), optimizer = Cbc.Optimizer)

##

optimize!(sp)

##

println("Objective value: $(objective_value(sp))")

println("PV: $(round(value.(sp[1, :u_pv]); digits = 2)) -- Wind: $(round(value.(sp[1, :u_wind]); digits = 2)) -- Battery: $(round(value.(sp[1, :u_storage]); digits = 2))")
println("Heat storage: $(round(value.(sp[1, :u_heat_storage]); digits = 2))")
fixed_investment = Dict((
    :u_pv => value.(sp[1, :u_pv]),
    :u_wind => value.(sp[1, :u_wind]),
    :u_storage => value.(sp[1, :u_storage]),
    :u_heatpump => value.(sp[1, :u_heatpump]),
    :u_heat_storage => value.(sp[1, :u_heat_storage])
))


# for s in 1:minimum((n, 10))
#     println("Second stage objective value in scenario $s: $(objective_value(sp, s))")
# end

##

plot_results(sp, pv, wind, demand, hd = heatdemand, s = 1, stage_1 = [:gci, :gco], stage_2 = [:gci2, :gco2])
##

n = 100
F_max = average_hourly_demand * 0.1 # Have ancillary services equal to 10% of our typical demand
t_max = length(pv) - es.parameters[2].defaults[:recovery_time]

scens = simple_flex_sampler(n, F_max, t_max)
##
es_fixed = define_energy_system(pv, wind, demand, heatdemand; p = pars, strict_flex = true)
sp_f = instantiate(es_fixed, scens, optimizer = Cbc.Optimizer)
fix_investment!(sp_f, string.(keys(fixed_investment)), values(fixed_investment))

optimize!(sp_f)
##

optimal_decision(sp_f)
plot_results(sp_f, pv, wind, demand, hd = heatdemand, s = 8, stage_1 = [:gci, :gco], stage_2 = [:gci2, :gco2], inv_dict = fixed_investment)
##
cost_pos, pot_pos, cost_neg, pot_neg = test_decision(sp, 1:t_max)

plot_flexibility(1:t_max, cost_pos, pot_pos, cost_neg, pot_neg, objective_value(sp))

##
cost_pos_f, pot_pos_f, cost_neg_f, pot_neg_f = test_decision(sp_f, 1:t_max)

plot_flexibility(1:t_max, cost_pos_f, pot_pos_f, cost_neg_f, pot_neg_f, objective_value(sp))

#label_av = ["+ flexibility in scenario-aware system", "+ flexibility in scenario-unaware system"]
begin
plt_av = plot();
flexibility_availability(plt_av, pot_pos, label = "positive flexbility unaware");
flexibility_availability(plt_av, pot_pos_f, label = "positive fleibility aware");
flexibility_availability(plt_av, pot_neg, label = "negative flexibility unaware");
flexibility_availability(plt_av, pot_neg_f, label = "negative flexibility aware");
display(plt_av)
end
##