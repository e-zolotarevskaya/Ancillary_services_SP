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
    costs = []
    decision = optimal_decision(p)
    for sign in [-1,1]
        for t in timesteps
            scenario = @scenario t_xi = t s_xi = sign F_xi = 150. probability = 1.
            c = evaluate_decision_wrapper(p, decision, scenario)
            append!(costs, c)
        end
    end
    N = length(timesteps)
    scenario_results = hcat(costs, vcat(timesteps, timesteps), vcat(-1*ones(N), ones(N)))
    return scenario_results
end

function test_decision_variate_F(p, timesteps; F_step = 50., F_max = 300.)
    decision = optimal_decision(p)
    costs = zeros(2*length(timesteps))
    F_potential = zeros(2*length(timesteps))
    for Sign in [-1,1]
        for t in timesteps
            i = t
            if Sign == 1
                i = t + length(timesteps)
            end
            for F in F_step:F_step:F_max
                scenario = @scenario t_xi = t s_xi = Sign F_xi = 150. probability = 1.
                c = evaluate_decision_wrapper(p, decision, scenario)
                if c!= Inf
                    F_potential[i] = F*Sign
                    costs[i] = c
                else
                    if costs[i] == 0
                        costs[i] = Inf
                        F_potential[i] = F*Sign
                        break
                    end
                end
            end
        end
    end
    print(F_potential)
    scenario_results = hcat(costs, vcat(timesteps, timesteps), F_potential)
    return scenario_results
end

function plot_flexibility(scen_results, ov)
    T = size(scen_results)[1]
    positive_flexibility = zeros(T ÷ 2)
    negative_flexibility = zeros(T ÷ 2)
    positive_potential = zeros(T ÷ 2)
    negative_potential = zeros(T ÷ 2)
    for i in 1:T
        if scen_results[i, 3] > 0
            if scen_results[i, 1] != Inf && scen_results[i, 1] != -Inf
                positive_flexibility[Int(scen_results[i, 2])] = scen_results[i, 1] - ov
                if scen_results != Inf
                    positive_potential[Int(scen_results[i, 2])] = scen_results[i, 3]
                end
            end
        else
            if scen_results[i, 1] != Inf && scen_results[i, 1] != -Inf
                negative_flexibility[Int(scen_results[i, 2])] = scen_results[i, 1] - ov
                if scen_results != -Inf
                    negative_potential[Int(scen_results[i, 2])] = scen_results[i, 3]
                end
            end
        end
    end
    plt_cost = plot()
    plt_pot = plot()
    plot!(plt_cost, positive_flexibility ./ positive_potential, label = "price of positive flexibility")
    plot!(plt_cost, negative_flexibility ./ negative_potential, label = "price of negative flexibility")
    plot!(plt_pot, positive_potential, fillrange = 0, fillalpha = 0.35, label = "positive flexibility potential")
    plot!(plt_pot, negative_potential, fillrange = 0, fillalpha = 0.35, label = "negative flexibility potential")
    display(plot(plt_cost, plt_pot, layout = (2, 1)))
    return positive_flexibility, negative_flexibility, positive_potential, negative_potential
end