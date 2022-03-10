using Plots

function plot_results(sp, pv, w, d; s=1, stage_1=[:gci, :gco], stage_2=[:gci2, :gco2], debug = false)
    plt_sto = plot()
    plt_invest = plot()
    plt = plot() # create separate plot for storage state of charge
    plot!(plt_invest, pv .* value(sp[1, :u_pv]), label="pv")
    plot!(plt_invest, w .* value(sp[1, :u_wind]), label="wind")
    plot!(plt_invest, d, label="demand")
    stor_flow_i = value.(sp[1, :sto_in])
    stor_flow_o = value.(sp[1, :sto_out])
    t_xi = scenarios(sp)[s].data.t_xi
    recovery_time = sp.stages[2].parameters[:recovery_time]
    stor_charge = [-sum(stor_flow_i[1:t])+sum(stor_flow_o[1:t]) for t in 1:length(pv)] .+ 0.5*value(sp[1, :u_storage])
    plot!(plt_sto, 1:length(pv), stor_charge, label="global storage charge")
    if debug
        #print("Maximum value of storage flow = "*string(maximum(value.(sp[1, :sto_in])))*"\n")
    end
    for var in stage_1
        plot!(plt, 1:length(pv), value.(sp[1, var]), label=string(var))
        if debug
            print("Maximum value of "*string(var)*" = "*string(maximum(value.(sp[1, var])))*"\n")
        end
    end
    #=for var in stage_2
        if debug
            print("Maximum value of "*string(var)*" = "*string(maximum(value.(sp[2, var], s)))*"\n")
        end
        plot!(plt, value.(sp[2, var], s).axes, value.(sp[2, var], s), label=string(var)*string(s), linestyle=:dash, linewidth=2)
    end=#
    for var in [:gci2, :gco2]
        plot!(plt, (t_xi+1):(t_xi+recovery_time), value.(sp[2, var], s), label=string(var)*string(s), linestyle=:dash, linewidth=2)
    end
    stor_flow_i2 = value.(sp[2, :sto_in2],s)
    stor_flow_o2 = value.(sp[2, :sto_out2],s)    
    stor_charge2 = [-sum(stor_flow_i2[1:t])+sum(stor_flow_o2[1:t]) for t in 1:(recovery_time+1)] .+ stor_charge[t_xi-1]
    plot!(plt_sto, (t_xi-1):(t_xi+recovery_time), vcat(stor_charge[t_xi-1],stor_charge2), label=string("stochastic storage charge")*string(s), linestyle=:dash, linewidth=2)
    if debug
        #print("Maximum value of sto2 = "*string(maximum(value.(sp[2, :sto2], s)))*"\n")
    end
    display(plot(plt_invest, plt, plt_sto, layout = (3,1)))
    #print(s)
end

function plot_flexibility(scen_results, ov)
    T = size(scen_results)[1]
    positive_flexibility = zeros(Int(T/2))
    negative_flexibility = zeros(Int(T/2))
    positive_potential = zeros(Int(T/2))
    negative_potential = zeros(Int(T/2))
    for i in 1:T
        if scen_results[i,3] > 0
            if scen_results[i,1] != Inf && scen_results[i,1] != -Inf
                positive_flexibility[Int(scen_results[i,2])] = scen_results[i,1] - ov
                if scen_results != Inf
                    positive_potential[Int(scen_results[i,2])] = scen_results[i, 3]
                end
            end
        else
            if scen_results[i,1] != Inf && scen_results[i,1] != -Inf
                negative_flexibility[Int(scen_results[i,2])] = scen_results[i,1] - ov
                if scen_results != -Inf
                    negative_potential[Int(scen_results[i,2])] = scen_results[i, 3]
                end
            end
        end
    end
    plt_cost = plot()
    plt_pot = plot()
    plot!(plt_cost, positive_flexibility./positive_potential, label = "price of positive flexibility")
    plot!(plt_cost, negative_flexibility./negative_potential, label = "price of negative flexibility")
    plot!(plt_pot, positive_potential, fillrange = 0, fillalpha = 0.35, label = "positive flexibility potential")
    plot!(plt_pot, negative_potential, fillrange = 0, fillalpha = 0.35, label = "negative flexibility potential")
    display(plot(plt_cost, plt_pot, layout = (2,1)))
    return positive_flexibility, negative_flexibility, positive_potential, negative_potential
end