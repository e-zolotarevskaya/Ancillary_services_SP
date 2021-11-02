cd(@__DIR__)
using Pkg
Pkg.activate(".")

using StochasticPrograms
using GLPK

##
@stochastic_model simple_model begin
    @stage 1 begin
        @decision(simple_model, x₁ >= 40)
        @decision(simple_model, x₂ >= 20)
        @objective(simple_model, Min, 100*x₁ + 150*x₂)
        @constraint(simple_model, x₁ + x₂ <= 120)
    end
    @stage 2 begin
        @uncertain q₁ q₂ d₁ d₂
        @recourse(simple_model, 0 <= y₁ <= d₁)
        @recourse(simple_model, 0 <= y₂ <= d₂)
        #@recourse(simple_model, [0,0].<=y.<=[d₁,d₂])
        @objective(simple_model, Max, q₁*y₁ + q₂*y₂)
        @constraint(simple_model, 6*y₁ + 10*y₂ <= 60*x₁)
        @constraint(simple_model, 8*y₁ + 5*y₂ <= 80*x₂)
        #= for i in 1:N 
            @constraint(simple_model, y[i]<a[i])
        end =#
        # @constraint(simple_model, [y[i]<a[i] for i in 1:N])
    end
end

ξ₁ = @scenario q₁ = 24.0 q₂ = 28.0 d₁ = 500.0 d₂ = 100.0 probability = 0.4

ξ₂ = @scenario q₁ = 28.0 q₂ = 32.0 d₁ = 300.0 d₂ = 300.0 probability = 0.6
##
sp = instantiate(simple_model, [ξ₁, ξ₂], optimizer = GLPK.Optimizer)

print(sp)

##
optimize!(sp)

objective_value(sp)

optimal_decision(sp)