cd(@__DIR__)
using Pkg
Pkg.activate(".")

using JuMP
using Cbc
using Plots
using CSV
using DataFrames

##


##Utility functions 
#=function clamped_sin(x)
    return sin(x-8)>0 ? 3. .*sin(x-8) : 0.
end=#
#pv = [clamped_sin(x/24. * 2. * pi) for x in timesteps]

#timesteps = 1:24
#pv = [clamped_sin(x/24. * 2. * pi) for x in timesteps]

#defining timesteps 
timesteps = 1:24

#pv, wind and demand data
demand = CSV.read("../data/demand_Industriepark.csv", DataFrame)[timesteps, 1]
wind = CSV.read("../data/wind_Karholz.csv", DataFrame)[timesteps, 1]
pv = CSV.read("../data/pv_Halle18.csv", DataFrame)[timesteps, 1]

#heatsystem demand assumption
heatdemand = zeros(Int, length(timesteps))
fill!(heatdemand, 100)

#Global system parameter 
c_i = .10 
c_o = .01
c_w = 1000.
c_pv = 1000.
c_storage = 400.
investment_capacity = 20000

c_heatstorage = 100


time_scale = 10. *365*24/length(timesteps)
storage_scale = 1


##midday without ancillary service 

#data supply
#demand = [1164.728, 1182.104, 1194.048]
#pv = [0.0391025641025641, 0.019230769230769232, 0.023717948717948717]
#wind = [0.9772296015180265, 0.9276462338090916, 0.9802821549377114]

m = Model(Cbc.Optimizer)

#varables declaration
@variable(m, 0 <= u_pv)
@variable(m, 0 <= u_w)
@variable(m, 0 <= u_storage)
@variable(m, gci[t in timesteps]>=0)
@variable(m, gco[t in timesteps]>=0)
@variable(m, storage[t in timesteps])
#variable heatsystem
#@variable(m, 0 <= u_heatpump)
@variable(m, 0 <= u_heatstorage)
@variable(m, heatstorage[t in timesteps])
@variable(m, heatpumpflow[t in timesteps])
#@variable(m, z[t in timesteps], Bin)

#energy balance
@constraint(m,[t in timesteps], gci[t]-gco[t]+u_pv*pv[t]+u_w*wind[t]-demand[t]+storage[t]-heatpumpflow[t]==0)
#without heat storage:  -heatdemand[t]* Coefficient of Performance (COP) from heatpump
#with heat storage: -(heatdemand[t]+heatstorage[t]) -> heatpumpflow[t]

#objective function
@objective(m, Min, c_pv/time_scale*u_pv + c_w/time_scale*u_w + u_storage*c_storage/time_scale + c_i*sum(gci) - c_o*sum(gco)) + + u_heatstorage*c_heatstorage/time_scale

#constraint for budget: 20.000
@constraint(m, u_pv + u_w + u_storage + u_heatstorage <= investment_capacity )


#defining battary constraint
#@constraint(m, [t in timesteps], -0.5*u_b*b_scale <= sum(b[1:t]))
#@constraint(m, [t in timesteps], sum(b[1:t]) <= 0.5*u_b*b_scale)
#@constraint(m, b[1]==b[length(timesteps)]) 

@constraint(m, [t in timesteps], -0.5*u_storage*storage_scale <= sum(storage[1:t]))
@constraint(m, [t in timesteps], sum(storage[1:t]) <= 0.5*u_storage*storage_scale)
@constraint(m, sum(storage) == 0)

#heatstorage constraints like in battary storage
@constraint(m, [t in timesteps], heatpumpflow[t] == heatdemand[t] + heatstorage[t])
@constraint(m, [t in timesteps], 0 <= sum(heatstorage[1:t])+ 0.5 * u_heatstorage)
@constraint(m, [t in timesteps], sum(heatstorage[1:t]) <= 0.5 * u_heatstorage)
@constraint(m, [t in timesteps], sum(heatstorage) == 0)

#defining heatpump constraint, with binary vaiable z
#@constraint(m, [t in timesteps], heatdemand[t] <= 1000 * (1 - z[t]) )
#@constraint(m, [t in timesteps], -heatpumpflow[t] * 1.1 + heatdemand[t]  == 1000*z[t])

print(m)

status = optimize!(m)

println("Objective value: ", getobjectivevalue(m))
println("u_pv = ", getvalue(u_pv))
println("u_w = ", getvalue(u_w))
println("u_storage = ", getvalue(u_storage))
println("u_heatstorage = ", getvalue(u_heatstorage))
println("gci =", getvalue(sum(gci))
##
