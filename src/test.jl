cd(@__DIR__)
using Pkg
Pkg.activate(".")

using StochasticPrograms
using GLPK
using Pipe
using Random
using Plots

## Utility functions
function random_time(Length;t=nothing)
    r_t = zeros(Length)
    xi = !isnothing(t) ? t : rand(1:Length)
    r_t[xi] = 1
    return r_t
end

function clamped_sin(x)
    return sin(x-8)>0 ? 3. .*sin(x-8) : 0.
end

function wind_timeseries(Length)
    return 4. .+ 2. .*rand(Length)
end
##
function demand_timeseries(Length)
    return 10. .+2. .*rand(Length)
end
function Theta(x)
    return x>=0 ? 1 : 0
end

# Global system parameter
c_i = 1
c_o = 1.2
n_y = 25
c_f = 1.3
F = 10
e_c = 5
timesteps = [i for i in 1:24]
pv = [clamped_sin(x/24. * 2. * pi) for x in timesteps]
w = wind_timeseries(24)
d = demand_timeseries(24)
##
function scenarios(n, timesteps)
    return [@scenario t_xi = random_time(length(timesteps)) s_xi = rand([1,-1]) probability = 1. /n]
end
##
#= Note on signs: s_xi<0 means that energy is requested, s_xi>0 means an additional consumption
To simplify the problem and avoid using inequalities in constraints, let's first study only s_xi<0
Then buy-back B>0, therefore it's associated with price c_i
=#
@stochastic_model em begin
    @stage 1 begin
        @decision(em, c_pv>=0)
        @decision(em, c_w>=0)
        @decision(em, c_b>=0)
        @decision(em, gci[t in timesteps]>=0)
        @decision(em, gco[t in timesteps]>=0)
        @objective(em, Min, c_pv+c_w+c_b+c_i*sum(gci)-c_o*sum(gco))
    end
    @stage 2 begin
        @uncertain t_xi s_xi #t_xi is a vector of zeros, with 1 at timestep of flexibility demand
        @recourse(em, b[t in timesteps])
        @recourse(em, B[t in timesteps])
        @constraint(em, [t in timesteps], b[t]+t_xi[t]*F*s_xi+B[t] == gco[t]-gci[t]-c_w*w[t]-c_pv*pv[t]+d[t]) # add buy-back function here
        @objective(em, Min, c_i*sum(B)) 
        @constraint(em, #=[t in timesteps],=# -0.5*c_b*e_c <= sum(b)) #initial charge = 0.5*c_b*e_c
        @constraint(em, sum(b)<=0.5*c_b*e_c) # not whole sum, only to t
    end
end 

## 

xi_1 = @scenario t_xi = random_time(24, t=3) s_xi = 1 probabibility = 0.5

xi_2 = @scenario t_xi = random_time(24, t=10) s_xi = -1 probabibility = 0.5
##
sp = instantiate(em, [xi_1, xi_2], optimizer = GLPK.Optimizer)

#@objective(model, Min, c'*x)

