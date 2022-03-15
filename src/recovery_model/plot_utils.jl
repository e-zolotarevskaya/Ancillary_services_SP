using Plots

function plot_results(sp, pv, w, d; s = 1, stage_1 = [:gci, :gco], stage_2 = [:gci2, :gco2])
    plt_sto = plot()
    plt_invest = plot()
    plt = plot() # create separate plot for storage state of charge
    plot!(plt_invest, pv .* value(sp[1, :u_pv]), label = "pv")
    plot!(plt_invest, w .* value(sp[1, :u_wind]), label = "wind")
    plot!(plt_invest, d, xlabel = "time, h", ylabel = "kWh", label = "demand")
    stor_flow_i = value.(sp[1, :sto_in])
    stor_flow_o = value.(sp[1, :sto_out])
    t_xi = scenarios(sp)[s].data.t_xi
    recovery_time = sp.stages[2].parameters[:recovery_time]
    stor_charge = [-sum(stor_flow_i[1:t]) + sum(stor_flow_o[1:t]) for t in 1:length(pv)] .+ 0.5 * value(sp[1, :u_storage])
    plot!(plt_sto, 1:length(pv), stor_charge, xlabel = "time, h", ylabel = "kWh", label = "background storage charge")

    for var in stage_1
        plot!(plt, 1:length(pv), value.(sp[1, var]), xlabel = "time, h", ylabel = "kWh", label = string(var))
    end

    for var in stage_2
        plot!(plt, (t_xi+1):(t_xi+recovery_time), value.(sp[2, var], s), xlabel = "time, h", ylabel = "kWh", label = string(var) * " in scenario " * string(s), linestyle = :dash, linewidth = 2)
    end
    stor_flow_i2 = value.(sp[2, :sto_in2], s)
    stor_flow_o2 = value.(sp[2, :sto_out2], s)
    stor_charge2 = [-sum(stor_flow_i2[1:t]) + sum(stor_flow_o2[1:t]) for t in 1:(recovery_time+1)] .+ stor_charge[t_xi-1]
    plot!(plt_sto, (t_xi-1):(t_xi+recovery_time), vcat(stor_charge[t_xi-1], stor_charge2), xlabel = "time, h", ylabel = "kWh", label = string("storage charge") * " in scenario " * string(s), linestyle = :dash, linewidth = 2)
    display(plot(plt_invest, plt, plt_sto, layout = (3, 1)))
    #print(s)
end

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
