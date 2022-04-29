using Plots
function evaluate_decision_wrapper(p, decision, scenario)
    cost = 0.
    try cost = evaluate_decision(p, decision, scenario)
    catch e
        cost = Inf
    end
    return cost
end

function test_decision(p, timesteps)
    decision = optimal_decision(p)
    L = length(timesteps)
    cost_pos_flex = zeros(L)
    cost_neg_flex = zeros(L)
    potential_pos_flex = zeros(L)
    potential_neg_flex = zeros(L)
    for t in timesteps
        potential_pos_flex[t], cost_pos_flex[t] = find_f_max(p,t,1,decision,10.)
        potential_neg_flex[t], cost_neg_flex[t] = find_f_max(p,t,-1,decision,10.)
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

function flexibility_availability(potential_pos, potential_neg)
    pot_pos = sort(unique(potential_pos))
    frac_pos = 1 .-[sum(potential_pos[:].<=f) for f in pot_pos]./length(potential_pos)
    pot_neg = sort(unique(potential_neg))
    frac_neg = [sum(potential_neg[:].<=f) for f in pot_neg]./length(potential_neg)

    #plt = plot(pot_neg, frac_neg)
    #plot!(plt, pot_pos, frac_pos)
    #display(plt)
    return frac_pos, frac_neg
end

function find_f_max(p,t,s,od,tol; maxiter = 100)
    a = 0.
    b = 10000.
    i = 0
    scen = @scenario t_xi = t s_xi = s F_xi = 0. probability = 1.
    cost = evaluate_decision_wrapper(p,od,scen)
    if cost == Inf
        return 0.,cost
    else
        while b-a>tol && i<maxiter
            i+=1
            scen = @scenario t_xi = t s_xi = s F_xi = (a+b)/2 probability = 1.
            cost = evaluate_decision_wrapper(p,od,scen)
            if cost == Inf
                b = (a+b)/2
            else
                local at = a
                a = (a+b)/2
                b = (3b-at)/2
            end
        end
        if cost == Inf
            scen = @scenario t_xi = t s_xi = s F_xi = (a+b)/2-tol probability = 1.
            cost = evaluate_decision_wrapper(p,od,scen)
        end
        println(i)
        return s*(a+b)/2,cost
    end
end

