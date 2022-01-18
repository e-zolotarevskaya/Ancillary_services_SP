# Gas problem von lecture 1, linear programming
cd(@__DIR__)
using Pkg
Pkg.activate(".")
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

#adding binary problem 

using Cbc


s = Model(Cbc.Optimizer) #achtung Cbc Solver! 


@variable(s, y, Bin) #Binäres Problem 
@variable(s, 0<= x_nl <= 28)
@variable(s, 0<= x_nor <= 27)
@variable(s, 0<= x_ru <= 35)
@variable(s, 0<=x_lng1 <= 15)
@variable(s, 0<= x_lng2 <= 12)

@objective(s, Min, 70x_nl + 104x_nor + 103x_ru + 105x_lng1 + 106x_lng2)
@constraint(s, x_nl + x_nor + x_ru + x_lng1 + x_lng2 >= 76)
@constraint(s, x_ru <= 1000*y) #Binäres Problem
@constraint(s, x_lng1 + x_lng2 <= 1000*(1 - y)) #Binäres Problem 
@constraint(s, x_nl + x_nor >= 0.5 * 76)

print(s)

status = optimize!(s)

println("Objective value: ", getobjectivevalue(s))