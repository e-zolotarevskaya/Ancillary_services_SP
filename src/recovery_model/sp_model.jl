using StochasticPrograms
using Random

##

default_es_pars = Dict((
    :c_i => .3,
    :c_o => .05,
    :c_sto_op => 0.0001,
    :asset_lifetime => 10.,
    :c_pv => 700.,
    :c_wind => 2000.,
    :c_storage => 600.,
    :inv_budget => 500000000.,
    :flexible_demand => 10.,
    :recovery_time => 72
))

function define_energy_system(pv, wind, demand; p = default_es_pars)
    number_of_hours = minimum([length(pv), length(demand), length(wind)])
    c_sto_op = p[:c_sto_op]
    c_i = p[:c_i]
    c_o = p[:c_o]
    recovery_time = p[:recovery_time]

    energy_system = @stochastic_model begin 
        @stage 1 begin
            @parameters begin
                # Expected lifetime of components, years
                asset_lifetime = p[:asset_lifetime]
                # Costs in Euro/1kWp
                c_pv = p[:c_pv]
                c_wind = p[:c_wind]
                c_storage = p[:c_storage]
                c_sto_op = c_sto_op
                c_i = c_i
                c_o = c_o
                # Euro
                inv_budget = p[:inv_budget] # Make the problem bounded
                # Shiftable demand, kW
                flexible_demand = p[:flexible_demand]
            end
            lifetime_factor = asset_lifetime * 365 * 24 / number_of_hours
            # Component units to be invested in, kWp
            @decision(model, u_pv >= 0)
            @decision(model, u_wind >= 0)
            @decision(model, u_storage >= 0)
            @constraint(model, c_pv * u_pv + c_wind * u_wind + c_storage * u_storage <= inv_budget)
            # Grid connection
            @decision(model, gci[t in 1:number_of_hours] >= 0)
            @decision(model, gco[t in 1:number_of_hours] >= 0)
            # Flexible demand
            @decision(model, fl_dem[t in 1:number_of_hours] >= 0)
            @constraint(model, sum(fl_dem) == flexible_demand)
            # Storage model
            @decision(model, sto_in[t in 1:number_of_hours] >= 0) # into the bus from storage
            @decision(model, sto_out[t in 1:number_of_hours] >= 0)
            @constraint(model, [t in 1:number_of_hours], -0.5 * u_storage <= -sum(sto_in[1:t]) + sum(sto_out[1:t]))
            @constraint(model, [t in 1:number_of_hours], -sum(sto_in[1:t]) + sum(sto_out[1:t]) <= 0.5 * u_storage)
            @constraint(model, -sum(sto_in) + sum(sto_out) == 0)
            # Energy balance
            @constraint(model, [t in 1:number_of_hours], 
            gci[t] - gco[t] + u_pv * pv[t] + u_wind * wind[t] - demand[t] - fl_dem[t] + sto_in[t] - sto_out[t] == 0)
            # Investment costs
            @objective(model, Min, (u_pv * c_pv + u_wind * c_wind + u_storage * c_storage) / lifetime_factor
            + c_i * sum(gci) - c_o * sum(gco) +
            c_sto_op * sum(sto_in) + c_sto_op * sum(sto_out))
        end
        @stage 2 begin
            @parameters begin
                recovery_time = recovery_time
                c_sto_op = c_sto_op
                c_i = c_i
                c_o = c_o
                
            end
            @uncertain t_xi s_xi F_xi # t_xi the time of flexibility demand, s_xi - sign (±1 or 0)
            @known(model, u_pv)
            @known(model, u_wind)
            @known(model, u_storage)
            @known(model, gci)
            @known(model, gco)
            @known(model, sto_in)
            @known(model, sto_out)
            @known(model, fl_dem)
            # Post event components
            @recourse(model, gci2[t in 1:recovery_time] >= 0)
            @recourse(model, gco2[t in 1:recovery_time] >= 0)
            @recourse(model, sto_in2[t in 1:(recovery_time + 1)] >= 0)
            @recourse(model, sto_out2[t in 1:(recovery_time + 1)] >= 0)
            @recourse(model, fl_dem2[t in 1:(recovery_time + 1)] >= 0)

            # Storage
            @constraint(model, [t in 1:(recovery_time + 1)], - 0.5 * u_storage <= - sum(sto_in[1:(t_xi - 1)]) + sum(sto_out[1:(t_xi - 1)]) - sum(sto_in2[1:t]) + sum(sto_out2[1:t]))
            @constraint(model, [t in 1:(recovery_time + 1)], - sum(sto_in[1:(t_xi - 1)]) + sum(sto_out[1:(t_xi - 1)]) - sum(sto_in2[1:t]) + sum(sto_out2[1:t]) <= 0.5 * u_storage)

            @constraint(model, -sum(sto_in2) + sum(sto_out2) == -sum(sto_in[t_xi:(t_xi + recovery_time)]) + sum(sto_out[t_xi:(t_xi + recovery_time)]))

            # Flexible demand
            @constraint(model, sum(fl_dem2) == sum(fl_dem[t_xi:t_xi + recovery_time]))

            # Event energy balance
            # The storage and other fast acting components use the recourse variables here.
            # They provide the balance. Grid connection is not allowed, as we are suporting the grid here. 
            @constraint(model, gci[t_xi] - gco[t_xi] + u_pv * pv[t_xi] + u_wind * wind[t_xi] - demand[t_xi] - fl_dem2[1] + sto_in2[1] - sto_out2[1] + F_xi * s_xi == 0)
            # Post event energy balance
            @constraint(model, [t in 1:recovery_time],
            gci2[t] - gco2[t]
            + u_pv * pv[t + t_xi] + u_wind * wind[t + t_xi]
            - demand[t + t_xi] - fl_dem2[t + 1] 
            + sto_in2[t + 1] - sto_out2[t + 1] == 0)
            @objective(model, Min,
            + c_i * (sum(gci2) - sum(gci[(t_xi + 1):(t_xi + recovery_time)]))
            - c_o * (sum(gco2) - sum(gco[(t_xi + 1):(t_xi + recovery_time)]))
            + c_sto_op * (sum(sto_in2) + sum(sto_out2) - sum(sto_in[t_xi:(t_xi + recovery_time)]) - sum(sto_out[t_xi:(t_xi + recovery_time)])))
        end
    end
    energy_system
end

function simple_flex_sampler(n, F_max, t_max)
    [@scenario t_xi = rand(1:t_max) s_xi = rand([-1, 1]) F_xi = rand() * F_max probability = 1/n 
        for i in 1:n]
end