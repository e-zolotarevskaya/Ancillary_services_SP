# Gas problem von lecture 1, linear programming
using Clp, JuMP
sources = ["NOR", "NL", "LNG1", "LNG2", "RU"]

# data from table
capacities = Dict(zip(sources, [27,28,15,12,35]))
production_cost = Dict(zip(sources, [54,65,88,88,36]))
transport_cost =Dict(zip(sources, [50,5,17,18,67]))
demand = 76 # 86 - 8 produced domestically

# model
m = Model(Clp.Optimizer)
@variable(m, x[sources] >= 0)
@constraint(m,[s in sources], x[s] <= capacities[s]) # capacity constraint
@constraint(m, sum(x[s] for s in sources) >= demand) # market clearing
@constraint(m, sum(x[s] for s in sources if s in ["NOR", "NL"]) >= demand*0.5) # 50 percent of demand met by EU or Norway
@objective(m, Min, sum(x[s]*(production_cost[s]+transport_cost[s]) for s in sources))
optimize!(m)
objective_value(m)