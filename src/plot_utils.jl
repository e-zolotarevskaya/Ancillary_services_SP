
function plot_results(sp, pv, w, d; s=1, stage_1=[:gci, :gco], stage_2=[:gci2, :gco2], debug = false)
    plt_sto = plot()
    plt_invest = plot()
    plt = plot() # create separate plot for storage state of charge
    plot!(plt_invest, pv .* value(sp[1, :u_pv]), label="pv")
    plot!(plt_invest, w .* value(sp[1, :u_wind]), label="wind")
    plot!(plt_invest, d, label="demand")
    stor_flow_i = value.(sp[1, :sto_in]).data
    stor_flow_o = value.(sp[1, :sto_out]).data

    stor_charge = [-sum(stor_flow_i[1:t])+sum(stor_flow_o[1:t]) for t in value.(sp[1, :sto_in]).axes[1]] .+ 0.5*value(sp[1, :u_storage])
    plot!(plt_sto, value.(sp[1, :sto_in]).axes, stor_charge, label="global storage charge")
    if debug
        #print("Maximum value of storage flow = "*string(maximum(value.(sp[1, :sto_in]).data))*"\n")
    end
    for var in stage_1
        plot!(plt, value.(sp[1, var]).axes, value.(sp[1, var]).data, label=string(var))
        if debug
            print("Maximum value of "*string(var)*" = "*string(maximum(value.(sp[1, var]).data))*"\n")
        end
    end
    for var in stage_2
        if debug
            print("Maximum value of "*string(var)*" = "*string(maximum(value.(sp[2, var], s).data))*"\n")
        end
        plot!(plt, value.(sp[2, var], s).axes, value.(sp[2, var], s).data, label=string(var)*string(s), linestyle=:dash, linewidth=2)
    end
    stor_flow_i = value.(sp[2, :sto_in2],s).data
    stor_flow_o = value.(sp[2, :sto_out2],s).data        
    stor_charge = [-sum(stor_flow_i[1:t])+sum(stor_flow_o[1:t]) for t in value.(sp[2, :sto_in2], s).axes[1]] .+ 0.5*value(sp[1, :u_storage])
    plot!(plt_sto, value.(sp[2, :sto_in2], s).axes, stor_charge, label=string("stochastic storage charge")*string(s), linestyle=:dash, linewidth=2)
    if debug
        #print("Maximum value of sto2 = "*string(maximum(value.(sp[2, :sto2], s).data))*"\n")
    end
    display(plot(plt_invest, plt, plt_sto, layout = (3,1)))
    #print(s)
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
en

function plot_difference(sp; s = 1)
    plt_gc = plot()
    plt_sto = plot()
    t_xi = scenarios(sp)[s].data.t_xi
    recovery_time = sp.stages[2].parameters[:recovery_time]

    plot!(plt_gc, -value.(sp[1, :gci])[t_xi+1:t_xi+recovery_time]+value.(sp[2, :gci2], s), xlabel = "time after the event, h", label = "gci2-gci")
    plot!(plt_gc, -value.(sp[1, :gco])[t_xi+1:t_xi+recovery_time]+value.(sp[2, :gco2], s), xlabel = "time after the event, h", label = "gco2-gco")
    
    stor_flow_i = value.(sp[1, :sto_in])[t_xi:t_xi+recovery_time]
    stor_flow_o = value.(sp[1, :sto_out])[t_xi:t_xi+recovery_time]
    stor_charge = [-sum(stor_flow_i[1:t]) + sum(stor_flow_o[1:t]) for t in 1:recovery_time+1] .+ 0.5 * value(sp[1, :u_storage])
    stor_flow_i2 = value.(sp[2, :sto_in2], s)
    stor_flow_o2 = value.(sp[2, :sto_out2], s)
    stor_charge2 = [-sum(stor_flow_i2[1:t]) + sum(stor_flow_o2[1:t]) for t in 1:(recovery_time+1)] .+ stor_charge[t_xi-1]
    
    plot!(plt_sto, stor_charge2-stor_charge, label = "soc2-soc")


    display(plot(plt_gc, plt_sto, layout = (2, 1)))
end