using Dates
function save_results(sp)
    timestamp = string(now())
    #filename = "optimization_results/results"*timestamp
    filename = "results"*timestamp

    s = scenarios(sp)
    n = length(s)
    n_recourse = length(optimal_recourse_decision(sp,1))
    recourse_results = zeros(n, n_recourse+2)
    for i in 1:n
        recourse_results[i,1] = s[i].data[:t_xi]
        recourse_results[i,2] = s[i].data[:s_xi]*s[i].data[:F_xi]
        recourse_results[i,3:end] = optimal_recourse_decision(sp,i)[:]
    end
    #=scens = zeros(n,2)
    scens[:,1] = [s[i].data[:t_xi] for i in 1:n]
    scens[:,2] = [s[i].data[:s_xi]*s[i].data[:F_xi] for i in 1:n]
    #scens[:,3] = [s[i].data[:t_xi]]
    df = DataFrame(scens, :auto)
    CSV.write(filename*"scenarios.csv")
    for i in 1:n
        ord = optimal_recourse_decision(sp, i)
        df = DataFrame(reshape(ord, (length(ord), 1)), :auto)
        CSV.write(filename*"scenario$i.csv", df |> DataFrame)
    end=#
    od = optimal_decision(sp)
    df = DataFrame(reshape(od, (length(od), 1)), :auto)
    CSV.write(filename*".csv", df |> DataFrame)
    df_recourse = DataFrame(reshape(recourse_results, (n, n_recourse+2)),:auto)
    CSV.write(filename*"recourse_results.csv", df_recourse |> DataFrame)
end
