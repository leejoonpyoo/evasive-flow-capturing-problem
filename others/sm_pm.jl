using JuMP, Gurobi

# Define and solve the SM model using JuMP
function solve_SM_model(g, A_f, N_f, F, distances, lamda_bars, c_f, w_e)
    f = length(F)
    n = nv(g)
    ε = 0.0001

    model = Model(Gurobi.Optimizer)

    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, u_f[1:f] >= 0)  # Non-negative variables
    @variable(model, π_k[1:f, 1:n] >= 0)  # Non-negative variables
    @variable(model, r_ij[1:f, 1:n, 1:n] >= 0)  # Non-negative variables

    # Objective function
    @objective(model, Min, sum(w_e * x[i, j] for i in 1:n, j in 1:n if has_edge(g, i, j)) +
                       sum(c_f * distances[i,j] * r_ij[a, i, j] for a in 1:f, i in 1:n, j in 1:n if (i, j) in A_f[F[a]]))

    # Constraints
    # for a in 1:f
    #     for b in 1:n
    #         if b ∉ N_f[F[a]]
    #             @constraint(model, π_k[a, b] == 0)
    #         end
    #     end
    # end

    for a in 1:f
        for i in 1:n
            for j in 1:n
                if (i, j) ∉ A_f[F[a]]
                    @constraint(model, r_ij[a, i, j] == 0)
                end
            end
        end
    end

    for i in 1:n
        for j in 1:n
            if !has_edge(g, i, j)
                @constraint(model, x[i, j] == 0)
            end
        end
    end

    # Constraints (25)
    for a in 1:f
        for (i, j) in A_f[F[a]]
            @constraint(model, π_k[a, j] <= π_k[a, i] + distances[i, j] + (lamda_bars[(F[a])] + ε)* x[i, j])
        end
    end

    # Constraints (26)
    for a in 1:f
        @constraint(model, (lamda_bars[(F[a])] + ε) * (1 - u_f[a]) + sum(distances[i, j] * r_ij[a, i, j] for (i, j) in A_f[F[a]]) <= π_k[a, F[a][2]] - π_k[a, F[a][1]])
    end

    # Constraints (27)
    for a in 1:f
        for i in N_f[F[a]]
            if i == F[a][1]  # Source node
                @constraint(model, sum(r_ij[a, i, j] for (i, j) in A_f[F[a]] if i == i) - sum(r_ij[a, j, i] for (j, i) in A_f[F[a]] if i == i) == u_f[a])
            elseif i == F[a][2]  # Destination node
                @constraint(model, sum(r_ij[a, i, j] for (i, j) in A_f[F[a]] if i == i) - sum(r_ij[a, j, i] for (j, i) in A_f[F[a]] if i == i) == -u_f[a])
            else  # Intermediate nodes
                @constraint(model, sum(r_ij[a, i, j] for (i, j) in A_f[F[a]] if i == i) - sum(r_ij[a, j, i] for (j, i) in A_f[F[a]] if i == i) == 0)
            end
        end
    end    

    # Constraints (28)
    for a in 1:f
        for (i, j) in A_f[F[a]]
            @constraint(model, r_ij[a, i, j] <= 1 - x[i, j])
        end
    end

    # Solve the model
    optimize!(model)

    return objective_value(model), value.(x), value.(u_f), value.(π_k), value.(r_ij)
end

# Define and solve the PM model using JuMP
function solve_PM_model(g, A_f, N_f, F, distances, lamda_bars, c_f, w_e)
    f = length(F)
    n = nv(g)

    model = Model(Gurobi.Optimizer)

    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, u_f[1:f] >= 0)  # Non-negative variables
    @variable(model, r_f[1:f, 1:n, 1:n] >= 0)  # Non-negative variables

    # Objective function
    @objective(model, Min, sum(w_e * x[i, j] for i in 1:n, j in 1:n if has_edge(g, i, j)) +
                       sum(c_f * distances[i,j] * r_f[a, i, j] for a in 1:f, i in 1:n, j in 1:n if (i, j) in A_f[F[a]]))

    # Constraints
    for a in 1:f
        for i in 1:n
            for j in 1:n
                if (i, j) ∉ A_f[F[a]]
                    @constraint(model, r_f[a, i, j] == 0)
                end
            end
        end
    end

    for i in 1:n
        for j in 1:n
            if !has_edge(g, i, j)
                @constraint(model, x[i, j] == 0)
            end
        end
    end

    # Constraints (37)
    for a in 1:f
        @constraint(model, 1 - sum(x[i, j] for (i, j) in A_f[F[a]]) <= u_f[a])
    end

    # Constraints (38)
    for a in 1:f
        for i in N_f[F[a]]
            if i == F[a][1]  # Source node
                @constraint(model, sum(r_f[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_f[a, j, i] for (j, i) in A_f[F[a]]) == u_f[a])
            elseif i == F[a][2]  # Destination node
                @constraint(model, sum(r_f[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_f[a, j, i] for (j, i) in A_f[F[a]]) == -u_f[a])
            else  # Intermediate nodes
                @constraint(model, sum(r_f[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_f[a, j, i] for (j, i) in A_f[F[a]]) == 0)
            end
        end
    end    

    # Constraints (39)
    for a in 1:f
        for (i, j) in A_f[F[a]]
            @constraint(model, r_f[a, i, j] <= 1 - x[i, j])
        end
    end

    # Constraints (40)
    for i in 1:n
        for j in 1:n
            @constraint(model, x[i, j] in Bin)
        end
    end

    # Constraints (41)
    for a in 1:f
        for (i, j) in A_f[F[a]]
            @constraint(model, u_f[a] * r_f[a, i, j] >= 0)
        end
    end

    # Solve the model
    optimize!(model)

    return objective_value(model), value.(x), value.(u_f), value.(r_f)
end