cd(@__DIR__)
using Pkg
Pkg.activate(".")

using StochasticPrograms
using GLPK
##
# Global system parameters
c_i = 1
c_o = 1.2
n_y = 25
c_f = 1.3
F = 10
e_c = 5
##

@stochastic_model em begin
    @stage 1 begin
        @decision(em, c_pv>=0)
        @decision(em, c_w>=0)
        @decision(em, c_b>=0)
        @decision(em, gci.>=0)
        @decision(em, gco.>=0)
        @objective(em, Min, c_pv+c_w+c_b+c_i*sum(gci)-c_o*sum(gco))
    end
    @stage 2 begin
        @uncertain t_xi s_xi #t_xi is a vector of zeros, with 1 at timestep of flexibility demand
        @recourse(em, gci-gco+c_w*w+c_pv*pv-d+b+F*s_xi*t_xi=0)
        @objective(em, Min, c_i*Theta(-B) - c_o*Theta(B))
        @constraint(em, sum(b)>=-0.5*c_b*e_c) #initial charge = 0.5*c_b*e_c
        @constraint(em, sum(b)<=0.5*c_b*e_c)
    end
end 



@objective(model, Min, c'*x)