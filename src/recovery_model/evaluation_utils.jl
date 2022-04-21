using Plots
function evaluate_decision_wrapper(p, decision, scenario)
    cost = 0.
    try cost = evaluate_decision(p, decision, scenario)
    catch e
        cost = Inf
    end
    return cost
end

function test_decision(p, timesteps; F_min = 100., F_step = 50., F_max = 300.)
    decision = optimal_decision(p)
    L = length(timesteps)
    cost_pos_flex = zeros(L)
    cost_neg_flex = zeros(L)
    potential_pos_flex = zeros(L)
    potential_neg_flex = zeros(L)
    for t in timesteps
        for F in F_min:F_step:F_max
            scenario = @scenario t_xi = t s_xi = 1 F_xi = F probability = 1.
            c = evaluate_decision_wrapper(p, decision, scenario)
            if c!= Inf
                potential_pos_flex[t] = F
                cost_pos_flex[t] = c
            else
                if cost_pos_flex[t] == 0
                    cost_pos_flex[t] = Inf
                    break
                end
            end
        end
        for F in F_min:F_step:F_max
            scenario = @scenario t_xi = t s_xi = -1 F_xi = F probability = 1.
            c = evaluate_decision_wrapper(p, decision, scenario)
            if c!= Inf
                potential_neg_flex[t] = -F
                cost_neg_flex[t] = c
            else
                if cost_neg_flex[t] == 0
                    cost_neg_flex[t] = Inf
                    break
                end
            end
        end
    end
    
    return cost_pos_flex, potential_pos_flex, cost_neg_flex, potential_neg_flex
end

function plot_flexibility(timesteps, cost_pos_flex, potential_pos_flex, cost_neg_flex, potential_neg_flex, obj_value)
    plt_cost = plot()
    plt_pot = plot()
    plot!(plt_cost, timesteps, (cost_pos_flex .- obj_value)./ potential_pos_flex, label = "price of positive flexibility")
    plot!(plt_cost, timesteps, (-cost_neg_flex .- obj_value)./ potential_neg_flex, label = "price of negative flexibility")
    plot!(plt_pot, timesteps, potential_pos_flex, fillrange = 0, fillalpha = 0.35, label = "positive flexibility potential")
    plot!(plt_pot, timesteps, potential_neg_flex, fillrange = 0, fillalpha = 0.35, label = "negative flexibility potential")
    display(plot(plt_cost, plt_pot, layout = (2, 1)))
end

function flexibility_availability(potential)
    pot = sort(unique(potential))
    #frac = [length(pot[pot.>=f]) for f in pot]./length(pot)
    frac = [sum(potential[:].<=f) for f in pot]./length(potential)
    plt = plot(pot, frac)
    display(plt)
end
