cd(@__DIR__)
using Pkg
Pkg.activate(".")

using StochasticPrograms
using GLPK
##

@stochastic_model em begin
    @stage 1 begin
        @decision(em,  c_pv>=0)
        @decision(em, c_w>=0)
        @decision(em, c_b>=0)
        @objective(em, Min, c_pv+c_w+c_b)
    end
    @stage 2 begin
        @uncertain t_xi s_xi
        @recourse(em, )
    end
end 