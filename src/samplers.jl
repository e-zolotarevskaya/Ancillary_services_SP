@define_scenario simple_scenario = begin
    t_xi::Int64
    s_xi::Int64
    @zero begin
        return simple_scenario(0, 0)
    end
    @expectation begin
        t_xi = Int(floor(sum([probability(s)*s.t_xi for s in scenarios])))
        s_xi = Int(floor(sum([probability(s)*s.s_xi for s in scenarios])))
        return simple_scenario(t_xi, s_xi)
    end
end

@sampler simple_sampler = begin
    timesteps::UnitRange{Int64}
    #n::Int64
    simple_sampler(timesteps::UnitRange{Int64}) = new(timesteps)
    
    @sample simple_scenario begin
        @parameters timesteps
        return simple_scenario(rand(timesteps[1:end-2]), rand([-1, 1]))
    end
end
## Try adding F to the scenario
@define_scenario F_scenario = begin
    t_xi::Int64
    s_xi::Int64
    F_xi::Float64
    @zero begin
        return F_scenario(0, 0, 0.)
    end
    @expectation begin
        t_xi = Int(floor(sum([probability(s)*s.t_xi for s in scenarios])))
        s_xi = Int(floor(sum([probability(s)*s.s_xi for s in scenarios])))
        F_xi = mean([probability(s)*s.F_xi for s in scenarios])
        return F_scenario(t_xi, s_xi, F_xi)
    end
end

@sampler F_sampler = begin
    timesteps::UnitRange{Int64}
    F_range::Tuple{Float64, Float64}
    #n::Int64
    F_sampler(timesteps::UnitRange{Int64}, F_range::Tuple{Float64, Float64}) = new(timesteps, F_range)
    
    @sample F_scenario begin
        @parameters timesteps F_range
        return F_scenario(rand(timesteps[1:end-2]), rand([-1, 1]), F_range[1]+rand()*F_range[2])
    end
end