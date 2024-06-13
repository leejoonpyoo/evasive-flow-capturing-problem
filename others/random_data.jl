using Graphs, Random, Statistics

function add_edge(g, u, v, max_degree)
    if !has_edge(g, u, v) && degree(g, u) < max_degree && degree(g, v) < max_degree
        add_edge!(g, u, v)
    end
end

function generate_25_node_network()
    Random.seed!(1234)

    # graph with 25 nodes and 43 edges
    g = SimpleGraph(25)

    # Ensure that each node has at least one edge
    for u in 1:25
        while degree(g, u) == 0
            v = rand(1:25)
            add_edge(g, u, v, 6)
        end
    end

    # Add remaining edges while maintaining the maximum degree constraint and total arc count
    while ne(g) < 43
        u, v = rand(1:25, 2)
        add_edge(g, u, v, 6)
    end

    # OD pairs
    # Initialize a matrix to store distances
    distances = fill(0, 25, 25)

    # Assign distances for existing edges
    for e in edges(g)
        u, v = src(e), dst(e)
        dist = rand(20:40)
        distances[u, v] = dist
        distances[v, u] = dist
    end

    # Assign distances for non-connected node pairs
    for u in 1:25
        for v in 1:25
            if u != v && distances[u, v] == 0
                dist = rand(40:380)
                distances[u, v] = dist
                distances[v, u] = dist
            end
        end
    end
    return g, distances
end