cd(@__DIR__)
using Pkg
Pkg.activate(".")

Pkg.add("Cbc")

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
timesteps = 1:3

#pv, wind and demand data
demand = CSV.read("../data/demand_Industriepark.csv", DataFrame)[timesteps, 1]
wind = CSV.read("../data/wind_Karholz.csv", DataFrame)[timesteps, 1]
pv = CSV.read("../data/pv_Halle18.csv", DataFrame)[timesteps, 1]

#Global system parameter 
c_i = .03 
c_o = .01
c_w = 1000.
c_pv = 1000.
c_b = 400.

time_scale = 10. *365*24/length(timesteps)
b_scale = .5


##midday without ancillary service 

#data supply
#demand = [1164.728, 1182.104, 1194.048]
#pv = [0.0391025641025641, 0.019230769230769232, 0.023717948717948717]
#wind = [0.9772296015180265, 0.9276462338090916, 0.9802821549377114]

m = Model(Cbc.Optimizer)

#varables declaration
@variable(m, 0 <= u_pv)
@variable(m, 0 <= u_w)
@variable(m, 0 <= u_b)
@variable(m, gci[t in timesteps]>=0)
@variable(m, gco[t in timesteps]>=0)
@variable(m, b[t in timesteps])

#energy balance
@constraint(m,[t in timesteps], gci[t]-gco[t]+u_pv*pv[t]+u_w*wind[t]-demand[t]+b[t]==0)

#objective function
@objective(m, Min, c_pv/time_scale*u_pv + c_w/time_scale*u_w + c_b/time_scale*u_b + c_i*sum(gci) - c_o*sum(gco))

#defining battary constraint
@constraint(m, [t in timesteps], -0.5*u_b*b_scale <= sum(b[1:t]))
@constraint(m, [t in timesteps], sum(b[1:t]) <= 0.5*u_b*b_scale)
@constraint(m, b[1]==b[length(timesteps)])


print(m)

status = optimize!(m)

println("Objective value: ", getobjectivevalue(m))
println("u_pv = ", getvalue(u_pv))
println("u_w = ", getvalue(u_w))

##
