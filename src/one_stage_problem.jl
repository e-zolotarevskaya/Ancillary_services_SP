cd(@__DIR__)
using Pkg
Pkg.activate(".")

using JuMP
using Cbc
using Plots
using CSV
using DataFrames
using GLPK

##


# Define the time interval and timesteps
t_s = 1
t_f = 168
offset = 0 #equals basic example offset = 6531 / pv and wind data 
timesteps = t_s:t_f

#pv, wind and demand data
#demand = CSV.read("../data/demand_Industriepark.csv", DataFrame)[timesteps, 1]
#wind = CSV.read("../data/wind_Karholz.csv", DataFrame)[timesteps, 1]
#pv = CSV.read("../data/pv_Halle18.csv", DataFrame)[timesteps, 1]
pv = CSV.read("../basic_example_normalized.csv", DataFrame)[timesteps.+offset, 3]
wind = CSV.read("../basic_example_normalized.csv", DataFrame)[timesteps.+offset, 4]
demand = CSV.read("../basic_example_normalized.csv", DataFrame)[timesteps.+offset, 2]


#heatsystem demand assumption
heatdemand = zeros(Int, length(timesteps))
fill!(heatdemand, 0)

#Global system parameter 
c_i = .30 
c_o = .05
c_w = 2000.
c_pv = 700.
c_storage = 600.
#investment_budget = 1000000.

#adding heat storage parameters 
c_heatstorage = 200.
c_heatpump = 10.
COP = 2. #coefficient of perfomance 


time_scale = 10. *365*24/length(timesteps)
storage_scale = 1



m = Model(GLPK.Optimizer)

#varables declaration
@variable(m, 0 <= u_pv)
@variable(m, 0 <= u_w)
@variable(m, 0 <= u_storage)
@variable(m, gci[t in timesteps]>=0)
@variable(m, gco[t in timesteps]>=0)
@variable(m, storage[t in timesteps])

#variable declaration heatsystem
@variable(m, 0 <= u_heatpump)
@variable(m, 0 <= u_heatstorage)
@variable(m, heatstorage[t in timesteps])
@variable(m, 0 <= heatpumpflow[t in timesteps])




#energy balance
@constraint(m,[t in timesteps], gci[t]-gco[t]+u_pv*pv[t]+u_w*wind[t]-demand[t]+storage[t] - heatpumpflow[t]/COP ==0) #add - heatpumpflow[t]/COP 

#objective function
@objective(m, Min, c_pv/time_scale*u_pv + c_w/time_scale*u_w + u_storage*c_storage/time_scale + c_i*sum(gci) - c_o*sum(gco) + u_heatstorage*c_heatstorage/time_scale + u_heatpump* c_heatpump/time_scale )  #add + u_heatstorage*c_heatstorage/time_scale + u_heatpump* c_heatpump/time_scale

#budget constraint
#@constraint(m, c_pv*u_pv + c_w*u_w + c_storage*u_storage + c_heatstorage*u_heatstorage + c_heatpump*u_heatpump  <= investment_budget )

#battery in the example
@constraint(m, [t in timesteps], -0.5*u_storage*storage_scale <= sum(storage[1:t]))
@constraint(m, [t in timesteps], sum(storage[1:t]) <= 0.5*u_storage*storage_scale)
@constraint(m, sum(storage) == 0)

#heat balance 
@constraint(m, [t in timesteps], heatpumpflow[t] + heatstorage[t] - heatdemand[t] == 0) #heatpumpflow can be either positive and negative
@constraint(m, [t in timesteps], heatpumpflow[t] <= u_heatpump) #necessary size of heatpump

#heatstorage constraints like in battary storage
@constraint(m, [t in timesteps], 0 <= sum(heatstorage[1:t])+ 0.5 * u_heatstorage)
@constraint(m, [t in timesteps], sum(heatstorage[1:t]) <= 0.5 * u_heatstorage)
@constraint(m, [t in timesteps], sum(heatstorage) == 0)




print(m)

status = optimize!(m)

println("Objective value: ", getobjectivevalue(m))
println("u_pv = ", getvalue(u_pv))
println("u_w = ", getvalue(u_w))
println("u_storage = ", getvalue(u_storage))
println("u_heatstorage = ", getvalue(u_heatstorage))
println("u_heatpump = ", getvalue(u_heatpump))
println("gci = ", value.(gci))
println("gco = ", value.(gco))




## Additional test code 

##midday without ancillary service 

#data supply
#demand = [1164.728, 1182.104, 1194.048]
#pv = [0.0391025641025641, 0.019230769230769232, 0.023717948717948717]
#wind = [0.9772296015180265, 0.9276462338090916, 0.9802821549377114]


##Utility functions 
#=function clamped_sin(x)
    return sin(x-8)>0 ? 3. .*sin(x-8) : 0.
end=#
#pv = [clamped_sin(x/24. * 2. * pi) for x in timesteps]

#timesteps = 1:24
#pv = [clamped_sin(x/24. * 2. * pi) for x in timesteps]

##
#defining battary constraint
#@constraint(m, [t in timesteps], -0.5*u_b*b_scale <= sum(b[1:t]))
#@constraint(m, [t in timesteps], sum(b[1:t]) <= 0.5*u_b*b_scale)
#@constraint(m, b[1]==b[length(timesteps)]) 

#defining heatpump constraint, with binary vaiable z
#@constraint(m, [t in timesteps], heatdemand[t] <= 1000 * (1 - z[t]) )
#@constraint(m, [t in timesteps], -heatpumpflow[t] * 1.1 + heatdemand[t]  == 1000*z[t])

#Energybalance without heat storage:  -heatdemand[t]* Coefficient of Performance (COP) from heatpump
#Energybalance with heat storage: -(heatdemand[t]+heatstorage[t]) -> heatpumpflow[t]