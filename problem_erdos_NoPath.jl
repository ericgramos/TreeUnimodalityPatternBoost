#Checks log concavity (location user input). Local search will run a fixed number of swaps (user input) before moving on. Skips optimization in the case that the score is already positive.

#Local optimization punishes the path graph.


include("constants.jl")
using JSON
using Polynomials
using DataStructures
using Random
const N = 56 #length of prufer code
#We will look for the breakage at alpha*n - beta
const alpha =  1/2
const beta  = 1
#maximum number of edge swaps allowed before optimization stops.
const max_swaps = 10
#determines whether the local optimization sorts edges before adding them.
const order_edges = 1

function star_list(n::Int)::Vector{Vector{Int}}
    adj = [Int[] for _ in 1:n]
    for i in 2:n
        push!(adj[1], i)
        push!(adj[i], 1)
    end
    return adj    
end

function prufer_to_adjacency_list(prufer::Vector{Int})::Vector{Vector{Int}}
    n = length(prufer) + 2
#if the machine generates a prufer code out of bounds, or generates a path, return a star graph with n vertices
    if any(v -> v ∉ 1:n, prufer) || length(Set(prufer)) == n-2
        return star_list(n)
    end
    
    degree = ones(Int, n)
    for i in prufer
        degree[i] += 1
    end

    # Use a min-heap of leaves
    heap = BinaryMinHeap{Int}()
    for v in 1:n
        if degree[v] == 1
            push!(heap, v)
        end
    end

    adjg = [Int[] for _ in 1:n]
    for v in prufer

        leaf = pop!(heap)
        push!(adjg[leaf], v)
        push!(adjg[v], leaf)
        degree[leaf] -= 1
        degree[v] -= 1

        if degree[v] == 1
            push!(heap, v)
        end
    end

    u = pop!(heap)
    v = pop!(heap)
    push!(adjg[u], v)
    push!(adjg[v], u)

    return adjg
end

function adjacency_list_to_prufer(adjg::Vector{Vector{Int}})::Vector{Int}
    n = length(adjg)
    # calculate the degree of each vertex by summing the rows
    degree = [length(adjg[i]) for i in 1:n]
    prufer = Int[]

    for x in 1:(n - 2)
        # Find the smallest leaf, and push its neighbor to prufer code
        leaf = findfirst(i -> degree[i] == 1, 1:n)
        if leaf === nothing
            # No leaf found, break out of the loop
            break
        end
        neighbor = first(adjg[leaf])
        push!(prufer, neighbor)

        # Remove the leaf
        filter!(x -> x != leaf, adjg[neighbor])
        degree[leaf] -= 1
        degree[neighbor] -= 1
        adjg[leaf] = Int[]  # Clear leaf's edges
    end

    return prufer
end

function find_cycle(adjg::Vector{Vector{Int}}, i::Int, j::Int)
    n = length(adjg)
    visited = falses(n)
    parent = fill(-1, n)

    # DFS to find path from i to j
    function dfs(u)
        visited[u] = true
        for v in adjg[u]
            if !visited[v]
                parent[v] = u
                if v == j || dfs(v)
                    return true
                end
            end
        end
        return false
    end

    dfs(i)

    # Reconstruct path from j to i
    path_nodes = Int[]
    curr = j
    while curr != -1
        push!(path_nodes, curr)
        curr = parent[curr]
    end
    reverse!(path_nodes)

    # Convert node path to edge list
    cycle_edges = [(path_nodes[k], path_nodes[k+1]) for k in 1:length(path_nodes)-1]
    push!(cycle_edges, (i, j)) 

    return cycle_edges
end

function IP(adjg::Vector{Vector{Int}})
    x = Polynomial([0, 1]) 
    I = Dict{Int, Polynomial}()
    V = Dict{Int, Polynomial}()
    N = Dict{Int, Polynomial}()
    C = [Int[] for _ in 1:length(adjg)]
    visited = Set{Int}()

    function find_inner_vertex(adjg::Vector{Vector{Int}})::Int
        for v in 1:length(adjg)
            if length(adjg[v]) >= 2
                return v
            end
        end
        return 1  # fallback
    end

    function IP_VISIT(adjg::Vector{Vector{Int}}, root::Int)
        push!(visited, root)
        if length(adjg[root]) == 1  
            V[root] = Polynomial([1])  
            N[root] = Polynomial([1])
            I[root] = Polynomial([1]) + x
            return
        end
        for u in adjg[root]
            if u ∉ visited
                push!(C[root], u)
                IP_VISIT(adjg, u)
            end
        end

        left = Polynomial([1])
        right = x
        for u in C[root]
            left *= I[u]
            right *= V[u]
        end
        V[root] = left
        N[root] = right
        I[root] = left + right
    end

    if length(adjg) == 1
        return Polynomial([1]) + x
    elseif length(adjg) == 2
        return Polynomial([1, 2])
    else
        root = find_inner_vertex(adjg)
        IP_VISIT(adjg, root)
        return I[root]
    end
end

"""prufer=[4,4,4,5]
adjg = prufer_to_adjacency_list(prufer)
poly = IP(adjg)
coeffs_array = coeffs(poly)
index = Int((N+2)/2) + 1
print(coeffs_array[index])"""

"""
- randomly add an edge (keep track of edges to add)
- find the cycle
- measure the log-concavity of the n/2 index of edges in the cycle
- find the edge that minimizes the difference
- check if it is the edge added, if yes, remove it and add another edge
- return the new graph
"""

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::OBJ_TYPE
    #input a database and a graph in prufer code, return the improved graph in prufer code

    # Check if current object has positive reward - if so, skip greedy search
    current_reward = reward_calc(obj)
    if current_reward > 0
        return obj  # Return the original object unchanged
    end

    prufer = obj
    adjg = prufer_to_adjacency_list(prufer)
    n = length(adjg)
    #degrees = map(length, adjg)
    
    none_edges = Set{Tuple{Int, Int}}()
    for i in 1:n
        for j in i+1:n
            if !(j in adjg[i])  
                push!(none_edges, (i, j))
            end
        end
    end

"""    sorted_edges = sort(
    collect(none_edges),
    by = kv -> (abs(kv[2][1] - kv[2][2]), -min(kv[2][1], kv[2][2])),
    rev = true)"""

"""    while !isempty(none_edges)
        # conitnue loop until we are able to decrease the goal
        # randomly add an edge
        edge_keys = collect(keys(none_edges))
        edge_to_add = edge_keys[rand(1:length(edge_keys))]
        # add an edge with maximum sum of degrees
        edge_to_add = argmax(none_edges)"""

    # Sort edges by having maximum degrees on both ends, if order_edges is set to true.
    if order_edges == 1
        none_edges = sort(collect(none_edges), by=kv -> (min(kv[2]...), max(kv[2]...)), rev=true)
    end

    swaps = 0

    for edge_to_add in none_edges
        #swaps will end our search as soon as enough optimizing moves have been made.
        if swaps == max_swaps
            break
        end

        i, j = edge_to_add
        push!(adjg[i], j)
        push!(adjg[j], i)

        # find the cycle created by the added edge
        cycle_edges = find_cycle(adjg, i, j)
        # loop through edges in the cycle and "measure" log-concavity at n/2 index
        cycle_log_count = Dict(edge => BigInt(0) for edge in cycle_edges)
        for (u,v) in cycle_edges        
            # Temporarily remove edge
            deleteat!(adjg[u], findfirst(==(v), adjg[u]))
            deleteat!(adjg[v], findfirst(==(u), adjg[v]))

            # Compute log-concavity at a(n/2-1) index
            poly = IP(adjg)
            coeffs_array = BigInt.(coeffs(poly))
            #TODO: checking at n/2-1 index
            index = floor(Int, alpha*(N+2)) - beta +1 
            
            # Pad coeffs_array with zeros if needed
            expected_length = index + 1
            if length(coeffs_array) < expected_length
                append!(coeffs_array, zeros(Int, expected_length - length(coeffs_array)))
            end
            
            a1 = coeffs_array[index-1]
            a2 = coeffs_array[index]
            a3 = coeffs_array[index+1]
            # update count
            cycle_log_count[(u,v)] = -(a2^2 - a1*a3)

            # Restore edge
            push!(adjg[u], v)
            push!(adjg[v], u)
        end

        #TODO
        # we select the edge that maximizes the score
        _, edge_to_delete = findmax(cycle_log_count)
        u, v = edge_to_delete
        filter!(x -> x != u, adjg[v])
        filter!(x -> x != v, adjg[u])

        #increment swaps, but decrement it if the edge you tried to delete is just the edge you added.
        swaps = swaps + 1

       if edge_to_delete == edge_to_add || edge_to_delete == (j, i)
            swaps = swaps - 1
       end
    end

    return adjacency_list_to_prufer(adjg)
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    prufer = obj
    adjg = prufer_to_adjacency_list(prufer)
    n = length(adjg)

    poly = IP(adjg)
    coeffs_array = BigInt.(coeffs(poly))
    # checking at floor(alpha*n) - beta index (sequence starts at 0)
    #TODO
    index = floor(Int, alpha*(N+2)) - beta +1
    
    # Pad coeffs_array with zeros if needed
    expected_length = index + 1
    if length(coeffs_array) < expected_length
        append!(coeffs_array, zeros(BigInt, expected_length - length(coeffs_array)))
    end

    a1 = coeffs_array[index-1]
    a2 = coeffs_array[index]
    a3 = coeffs_array[index+1]

    degrees = map(length, adjg)
    if count(==(1), degrees) == a1 - 1 && count(==(a1 - 1), degrees) == 1
        return -1e9  # punishment for star graph
    end  

    return -(a2^2 - a1*a3)
end

function prufer_to_json_string(prufer::Vector{Int})::String
    return JSON.json(prufer)
end

function json_string_to_prufer(json_str::String)::Vector{Int}
    return JSON.parse(json_str)
end

function count_leaves(prufer::Vector{Int})
    n = length(prufer) + 2
    prufer_set = Set(prufer)
    leaves = [v for v in 1:n if v ∉ prufer_set]
    return length(leaves)
end

function empty_starting_point()::OBJ_TYPE
    #TODO: recover sample generation   
    prufer = rand(1:(N+2), N)
    # check if the number of leaves is less than or equal to n/2
    #if count_leaves(prufer) <= 5
    """if count_leaves(prufer) >= floor(Int, (N+2)/2) || count_leaves(prufer) <= 5
        return empty_starting_point()
    # if not, generate a new prufer code
    else
        return prufer
    end"""
    return prufer
end
