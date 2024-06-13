using Graphs, Random, DataStructures

################################## Data Generation ##################################
function add_edge(g, u, v, distance, distances)
    if !has_edge(g, u, v)
        add_edge!(g, u, v)
        distances[u, v] = distance
        distances[v, u] = distance
    end
end

function generate_25_node_network()
    g = SimpleGraph(25)

    # Initialize a matrix to store lengths and distances
    distances = fill(0, 25, 25)

    # Add edges and lengths according to the provided graph
    edges_with_distances = [
        (1, 2, 40), (1, 3, 30), (1, 4, 40), (1, 5, 60),
        (2, 6, 30), (2, 7, 50), (2, 8, 70),
        (3, 9, 60), (3, 10, 50), (3, 11, 40),
        (4, 12, 40), (4, 13, 40),
        (5, 14, 70), (5, 15, 80),
        (6, 16, 60), (6, 17, 70),
        (7, 18, 80), (7, 19, 90),
        (8, 20, 90), (8, 21, 70),
        (9, 22, 80), (9, 23, 70),
        (10, 24, 90), (10, 25, 40),
        (11, 22, 30), (12, 23, 30),
        (13, 24, 30), (14, 25, 30)
    ]

    for (u, v, distance) in edges_with_distances
        add_edge(g, u, v, distance, distances)
    end

    return g, distances
end

##################################### Functions #####################################
function shortest_path_distances(g, start, goal, distances)
    if has_edge(g, start, goal)
        return distances[start, goal]
    else
        println("Calculating shortest path length from ", start, " to ", goal)
        dijkstra_distances = dijkstra_shortest_paths(g, start, distances)
        return dijkstra_distances[goal]
    end
end

function dijkstra_shortest_paths(g::Graph, start::Int, distances::Matrix{Int})
    n = nv(g)
    dist = fill(typemax(Int), n)
    dist[start] = 0
    pq = PriorityQueue{Int, Int}()
    enqueue!(pq, start => 0)

    while !isempty(pq)
        current_node, current_dist = peek(pq)
        dequeue!(pq)

        if current_dist > dist[current_node]
            continue
        end

        for neighbor in neighbors(g, current_node)
            edge_weight = distances[current_node, neighbor]
            if edge_weight > 0
                new_dist = dist[current_node] + edge_weight
                if new_dist < dist[neighbor]
                    dist[neighbor] = new_dist
                    if haskey(pq, neighbor)
                        delete!(pq, neighbor)
                    end
                    enqueue!(pq, neighbor => new_dist)
                end
            end
        end
    end

    return dist
end

function find_paths_within_tolerance(g, distances, lamda, start, goal)
    paths = []
    stack = [([start], 0)]  # Stack holds tuples of (path, total_length)
    base_distance = shortest_path_distances(g, start, goal, distances)
    tolerance = lamda * base_distance
    println("tolerance: ", tolerance)

    while !isempty(stack)
        path, current_distance = pop!(stack)
        current_node = path[end]

        if current_node == goal
            push!(paths, path)
        else
            for neighbor in neighbors(g, current_node)
                if neighbor ∉ path
                    new_path = vcat(path, [neighbor])
                    new_dist = current_distance + distances[current_node, neighbor]
                    if new_dist <= tolerance
                        push!(stack, (new_path, new_dist))
                    end
                end
            end
        end
    end

    return paths
end

function calculate_path_cost(path, distances, c_f)
    return sum(c_f * distances[path[i], path[i+1]] for i in 1:(length(path)-1))
end

function get_arc_set(path)
    return [(path[i], path[i+1]) for i in 1:(length(path)-1)]
end

function find_paths_and_calculate_costs(g, distances, lamda, c_f)
    F = []
    P = Dict{Tuple{Int, Int}, Vector{Vector{Int}}}()
    path_costs = Dict{Tuple{Int, Int}, Vector{Float64}}()
    A = Dict{Tuple{Int, Int}, Vector{Vector{Tuple{Int, Int}}}}()

    for u in 1:25
        for v in (u+1):25
            push!(F, (u, v))
            paths = find_paths_within_tolerance(g, distances, lamda, u, v)
            println(paths)
            P[(u, v)] = paths
            A[(u, v)] = [get_arc_set(path) for path in paths]
            path_costs[(u, v)] = [calculate_path_cost(path, distances, c_f) for path in paths]
        end
    end

    return F, P, path_costs, A
end

function lamda_bar(g, distances, lamda)
    F = []
    lamda_bars = Dict{Tuple{Int, Int}, Float64}()
    for u in 1:25
        for v in (u+1):25
            push!(F, (u, v))
            shortest_distance = dijkstra_shortest_paths(g, u, distances)[v]
            lamda_bars[(u, v)] = lamda * shortest_distance
        end
    end
    return F, lamda_bars 
end

function calculate_nondominant_arcs(g, distances, lamda)
    F, lamda_bars = lamda_bar(g, distances, lamda)
    A_f = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()

    for (u, v) in F
        lambda_bar_f = lamda_bars[(u, v)]
        A_f[(u, v)] = []

        for edge in edges(g)
            i, j = src(edge), dst(edge)

            ξ_s_f_i = dijkstra_shortest_paths(g, u, distances)[i]
            ξ_j_t_f = dijkstra_shortest_paths(g, j, distances)[v]
            d_ij = distances[i, j]

            if ξ_s_f_i + d_ij + ξ_j_t_f <= lambda_bar_f
                push!(A_f[(u, v)], (i, j))
            end
        end
    end
    return A_f
end

function calculate_nondominant_nodes(A_f)
    N_f = Dict{Tuple{Int, Int}, Set{Int}}()
    for (u, v) in keys(A_f)
        N_f[(u, v)] = Set{Int}()
        for (i, j) in A_f[(u, v)]
            push!(N_f[(u, v)], i)
            push!(N_f[(u, v)], j)
        end
    end
    return N_f
end

# function shortest_path_dijkstra(g::Graph, start::Int, goal::Int, distances::Matrix{Int})
#     n = nv(g)
#     dist = fill(typemax(Int), n)
#     dist[start] = 0
#     prev = fill(-1, n)  # To reconstruct the path
#     pq = PriorityQueue{Int, Int}()  # node => distance
#     enqueue!(pq, start => 0)

#     while !isempty(pq)
#         current_node, current_dist = peek(pq)
#         dequeue!(pq)

#         if current_node == goal
#             break
#         end

#         for neighbor in neighbors(g, current_node)
#             edge_weight = distances[current_node, neighbor]
#             if edge_weight > 0  # Ignore zero-length edges (non-edges)
#                 new_dist = current_dist + edge_weight
#                 if new_dist < dist[neighbor]
#                     dist[neighbor] = new_dist
#                     prev[neighbor] = current_node
#                     enqueue!(pq, neighbor => new_dist)
#                 end
#             end
#         end
#     end

#     # Reconstruct the shortest path
#     path = []
#     u = goal
#     while u != -1
#         push!(path, u)
#         u = prev[u]
#     end
#     reverse!(path)  # The path is reconstructed backwards

#     return path
# end

function shortest_path_dijkstra(g::Graph, start::Int, goal::Int, distances::Matrix{Int})
    n = nv(g)
    dist = fill(typemax(Int), n)
    dist[start] = 0
    prev = fill(-1, n)  # To reconstruct the path
    pq = PriorityQueue{Int, Int}()  # node => distance
    enqueue!(pq, start => 0)

    while !isempty(pq)
        current_node, current_dist = peek(pq)
        dequeue!(pq)

        if current_node == goal
            break
        end

        if current_dist > dist[current_node]
            continue
        end

        for neighbor in neighbors(g, current_node)
            edge_weight = distances[current_node, neighbor]
            if edge_weight > 0  # Ignore zero-length edges (non-edges)
                new_dist = dist[current_node] + edge_weight
                if new_dist < dist[neighbor]
                    dist[neighbor] = new_dist
                    prev[neighbor] = current_node
                    if haskey(pq, neighbor)
                        delete!(pq, neighbor)
                    end
                    enqueue!(pq, neighbor => new_dist)
                end
            end
        end
    end

    # Reconstruct the shortest path
    path = []
    u = goal
    while u != -1
        push!(path, u)
        u = prev[u]
    end
    reverse!(path)  # The path is reconstructed backwards

    return path
end

function find_shortest_path_arcs(g, F, distances)
    shortest_path_arcs = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()
    for (u, v) in F
        path = shortest_path_dijkstra(g, u, v, distances)
        arcs = [(path[i], path[i+1]) for i in 1:(length(path)-1)]
        shortest_path_arcs[(u, v)] = arcs
    end
    return shortest_path_arcs
end