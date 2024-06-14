include("functions.jl")
include("model.jl")

using Dates  

################################## 0. Set data ##################################
g, distances = generate_25_node_network()

# input parameters
c_f = 0.2  # [0.025, 0.05, 0.10, 0.20]                      # damage costs (dollars per kilometer)
w_e = 360  # [10000, 60000, 110000, 160000, 260000, 360000]    # opening costs (thousands of dollars)
lamda = 2.0  # [1.0, 1.5, 2.0]                                # Deviation tolerance

################################## 1. MRS model ##################################
F, P, path_costs, A = find_paths_and_calculate_costs(g, distances, lamda, c_f)

start_time = now()
obj_mrs, x_mrs, y_f_mrs, y_p_mrs, z_mrs = MRS(g, F, P, path_costs, A, w_e)
end_time = now()

println("MRS Model:")
println("Objective value: ", obj_mrs)
print_time_difference(start_time, end_time)

################################## 2. SM model ##################################
F, lamda_bars = lamda_bar(g, distances, lamda)
A_f = calculate_nondominant_arcs(g, distances, lamda)
N_f = calculate_nondominant_nodes(A_f)

start_time = now()
obj_sm, x_sm, u_f_sm, π_k_sm, r_ij_sm = SM(g, A_f, N_f, F, distances, lamda_bars, c_f, w_e)
end_time = now()

println("SM Model:")
println("Objective value: ", obj_sm)
print_time_difference(start_time, end_time)

################################## 3. PM model ##################################
start_time = now()
model, obj_pm, x_pm, u_f_pm, r_fij_pm = PM(g, A_f, N_f, F, distances, c_f, w_e, P, A)
end_time = now()

println("PM Model:")
println("Objective value: ", obj_pm)
print_time_difference(start_time, end_time)

################################## 4. PM cutting plane model ##################################
shortest_path_arcs = find_shortest_path_arcs(g, F, distances)

start_time = now()
model, obj_pm_cp, x_pm_cp, u_f_pm_cp, r_fij_pm_cp = PM_cp(g, A_f, N_f, F, distances, c_f, w_e, P, A, shortest_path_arcs)
end_time = now()

println("PM Cutting Plane Model:")
println("Objective value: ", obj_pm_cp)
print_time_difference(start_time, end_time)