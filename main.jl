include("functions.jl")
include("model.jl")

################################## 0. Set data ##################################
g, distances = generate_25_node_network()

# input parameters
c_f = 0.025  # [0.025, 0.05, 0.10, 0.20]                      # damage costs (dollars per kilometer)
w_e = 10  # [10000, 60000, 110000, 160000, 260000, 360000]    # opening costs (thousands of dollars)
lamda = 1.0  # [1.0, 1.5, 2.0]                                # Deviation tolerance

################################## 1. MRS model ##################################
F, P, path_costs, A = find_paths_and_calculate_costs(g, distances, lamda, c_f)
# compile gurobi
obj_mrs, x_mrs, y_f_mrs, y_p_mrs, z_mrs = MRS(g, F, P, path_costs, A, w_e)
obj_mrs, x_mrs, y_f_mrs, y_p_mrs, z_mrs = MRS(g, F, P, path_costs, A, w_e)

println("Objective value: ", obj_mrs)
println("x: ", x_mrs)
println("y_f: ", y_f_mrs)
println("y_p: ", y_p_mrs)
println("z: ", z_mrs)

################################## 2. SM model ##################################
F, lamda_bars = lamda_bar(g, distances, lamda)
A_f = calculate_nondominant_arcs(g, distances, lamda)
N_f = calculate_nondominant_nodes(A_f)
obj_sm, x_sm, u_f_sm, π_k_sm, r_ij_sm = SM(g, A_f, N_f, F, distances, lamda_bars, c_f, w_e)

println("Objective value: ", obj_sm)
# println("x: ", x_sm)
# println("u_f: ", u_f_sm)
# println("π_k: ", π_k_sm)
# println("r_ij: ", r_ij_sm)

################################## 3. PM model ##################################
model, obj_pm, x_pm, u_f_pm, r_fij_pm = PM(g, A_f, N_f, F, distances, c_f, w_e, P, A)

################################## 4. PM cutting plane model ##################################
shortest_path_arcs = find_shortest_path_arcs(g, F, distances)
model, obj_pm, x_pm, u_f_pm, r_fij_pm = PM_cp(g, A_f, N_f, F, distances, c_f, w_e, P, A, shortest_path_arcs)








