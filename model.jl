include("functions.jl")

using JuMP, Gurobi, CPLEX

function MRS(g, F, P, path_costs, A, w_e)
    n = nv(g)
    f = length(F)
    max_p = maximum(length(P[F[i]]) for i in 1:f)

    # Initialize the optimization model
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)

    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, y_f[1:f], Bin)
    @variable(model, y_p[1:f, 1:max_p], Bin)
    @variable(model, z[1:f, 1:max_p], Bin)

    # Objective function
    @objective(model, Min, sum(w_e * x[i, j] for i in 1:n, j in 1:n if has_edge(g, i, j)) + sum(path_costs[F[a]][b] * z[a, b] for a in 1:f, b in 1:max_p if b <= length(P[F[a]])))

    # Basic constraints
    for a in 1:f
        for b in 1:max_p
            if b > length(P[F[a]])
                @constraint(model, y_p[a, b] == 0)
                @constraint(model, z[a, b] == 0)
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

    # Additional constraints
    for a in 1:f
        for b in 1:length(P[F[a]])
            @constraint(model, y_p[a, b] <= sum(x[i,j] for (i,j) in A[F[a]][b]))  # (2)
            for (i, j) in A[F[a]][b]
                @constraint(model, x[i, j] <= y_p[a, b])  # (3)
            end
            @constraint(model, y_f[a] <=  y_p[a, b] )  # (4)
            @constraint(model, z[a, b] <= 1 - y_p[a, b])  # (5)
        end
        @constraint(model, 1 - y_f[a] <= sum(z[a, b] for b in 1:length(P[F[a]])))  # (6)
    end

    # Solve the model
    optimize!(model)

    return objective_value(model), value.(x), value.(y_f), value.(y_p), value.(z)
end

function SM(g, A_f, N_f, F, distances, lamda_bars, c_f, w_e)
    f = length(F)
    n = nv(g)
    ε = 0.0001

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)

    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, u_f[1:f] >= 0)  # Non-negative variables
    @variable(model, π_k[1:f, 1:n] >= 0)  # Non-negative variables
    @variable(model, r_fij[1:f, 1:n, 1:n] >= 0)  # Non-negative variables

    # Objective function
    @objective(model, Min, sum(w_e * x[i, j] for i in 1:n, j in 1:n if has_edge(g, i, j)) +
                       sum(c_f * distances[i,j] * r_fij[a, i, j] for a in 1:f, i in 1:n, j in 1:n if (i, j) in A_f[F[a]]))

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
                    @constraint(model, r_fij[a, i, j] == 0)
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
        @constraint(model, (lamda_bars[(F[a])] + ε) * (1 - u_f[a]) + sum(distances[i, j] * r_fij[a, i, j] for (i, j) in A_f[F[a]]) <= π_k[a, F[a][2]] - π_k[a, F[a][1]])
    end

    # Constraints (27)
    for a in 1:f
        for i in N_f[F[a]]
            if i == F[a][1]  # Source node
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]] if i == i) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]] if i == i) == u_f[a])
            elseif i == F[a][2]  # Destination node
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]] if i == i) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]] if i == i) == -u_f[a])
            else  # Intermediate nodes
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]] if i == i) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]] if i == i) == 0)
            end
        end
    end    

    # Constraints (28)
    for a in 1:f
        for (i, j) in A_f[F[a]]
            @constraint(model, r_fij[a, i, j] <= 1 - x[i, j])
        end
    end

    # Solve the model
    optimize!(model)

    return objective_value(model), value.(x), value.(u_f), value.(π_k), value.(r_fij)
end

function PM(g, A_f, N_f, F, distances, c_f, w_e, P, A)
    f = length(F)
    n = nv(g)

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)

    @variable(model, x[1:n, 1:n], Bin)          # Constraints (40)
    @variable(model, u_f[1:f] >= 0)             # Constraints (41)
    @variable(model, r_fij[1:f, 1:n, 1:n] >= 0)   # Constraints (41)

    # Objective function
    @objective(model, Min, sum(w_e * x[i, j] for i in 1:n, j in 1:n if has_edge(g, i, j)) +
                       sum(c_f * distances[i,j] * r_fij[a, i, j] for a in 1:f, i in 1:n, j in 1:n if (i, j) in A_f[F[a]]))

    # Constraints
    for a in 1:f
        for i in 1:n
            for j in 1:n
                if (i, j) ∉ A_f[F[a]]
                    @constraint(model, r_fij[a, i, j] == 0)
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
        for p in 1:length(P[F[a]])
            arc_set = A[F[a]][p]
            @constraint(model, 1 - sum(x[i, j] for (i, j) in arc_set) <= u_f[a])
        end
    end

    # Constraints (38)
    for a in 1:f
        for i in N_f[F[a]]
            if i == F[a][1]  # Source node
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]]) == u_f[a])
            elseif i == F[a][2]  # Destination node
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]]) == -u_f[a])
            else  # Intermediate nodes
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]]) == 0)
            end
        end
    end    

    # Constraints (39)
    for a in 1:f
        for (i, j) in A_f[F[a]]
            @constraint(model, r_fij[a, i, j] <= 1 - x[i, j])
        end
    end

    # Solve the model
    optimize!(model)

    return model, objective_value(model), value.(x), value.(u_f), value.(r_fij)
end

function PM_cp(g, A_f, N_f, F, distances, c_f, w_e, P, A, shortest_path_arcs)
    f = length(F)
    n = nv(g)

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)

    @variable(model, x[1:n, 1:n], Bin)          # Constraints (40)
    @variable(model, u_f[1:f] >= 0)             # Constraints (41)
    @variable(model, r_fij[1:f, 1:n, 1:n] >= 0)   # Constraints (41)

    # Objective function
    @objective(model, Min, sum(w_e * x[i, j] for i in 1:n, j in 1:n if has_edge(g, i, j)) +
                       sum(c_f * distances[i,j] * r_fij[a, i, j] for a in 1:f, i in 1:n, j in 1:n if (i, j) in A_f[F[a]]))

    # Constraints
    for a in 1:f
        for i in 1:n
            for j in 1:n
                if (i, j) ∉ A_f[F[a]]
                    @constraint(model, r_fij[a, i, j] == 0)
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

    # Constraints (37) initially considering only the shortest path
    for a in 1:f
        spa = shortest_path_arcs[F[a]]
        @constraint(model, 1 - sum(x[i, j] for (i, j) in spa) <= u_f[a])
    end

    # Constraints (38)
    for a in 1:f
        for i in N_f[F[a]]
            if i == F[a][1]  # Source node
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]]) == u_f[a])
            elseif i == F[a][2]  # Destination node
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]]) == -u_f[a])
            else  # Intermediate nodes
                @constraint(model, sum(r_fij[a, i, j] for (i, j) in A_f[F[a]]) - sum(r_fij[a, j, i] for (j, i) in A_f[F[a]]) == 0)
            end
        end
    end    

    # Constraints (39)
    for a in 1:f
        for (i, j) in A_f[F[a]]
            @constraint(model, r_fij[a, i, j] <= 1 - x[i, j])
        end
    end

    function add_cutting_planes!(model, x_val, u_f_val, F, A_f, P, A)
        added_constraints = false
        for a in 1:length(F)
            for p in 1:length(P[F[a]])
                arc_set = A[F[a]][p]
                if sum(x_val[i, j] for (i, j) in arc_set) < 1 - u_f_val[a] - 1e-6
                    @constraint(model, 1 - sum(x[i, j] for (i, j) in arc_set) <= u_f[a])
                    added_constraints = true
                end
            end
        end
        return added_constraints
    end

    optimize!(model)
    
    while true
        x_val = value.(x)
        u_f_val = value.(u_f)
        if !add_cutting_planes!(model, x_val, u_f_val, F, A_f, P, A)
            break
        end
        optimize!(model)
    end

    return model, objective_value(model), value.(x), value.(u_f), value.(r_fij)
end