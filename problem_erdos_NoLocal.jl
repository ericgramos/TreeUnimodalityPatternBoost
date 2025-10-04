#This version of the code does not do any local optimization. The greedy search function just returns the input unchanged.


include("constants.jl")
using JSON
using Polynomials
using DataStructures
using Random
const N = 58 #length of prufer code
#We will look for the breakage at alpha*(N+2) - beta
const alpha =  1/2
const beta  = 0



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
#if the machine generates a prufer code out of bounds, return a star graph with n vertices
    if any(v -> v ∉ 1:n, prufer)
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
    x = Polynomial([BigInt(0), BigInt(1)]) 
    I = Dict{BigInt, Polynomial}()
    V = Dict{BigInt, Polynomial}()
    N = Dict{BigInt, Polynomial}()
    C = [BigInt[] for _ in 1:length(adjg)]
    visited = Set{BigInt}()

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
            V[root] = Polynomial([BigInt(1)])  
            N[root] = Polynomial([BigInt(1)])
            I[root] = Polynomial([BigInt(1)]) + x
            return
        end
        for u in adjg[root]
            if u ∉ visited
                push!(C[root], u)
                IP_VISIT(adjg, u)
            end
        end

        left = Polynomial([BigInt(1)])
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
        return Polynomial([BigInt(1)]) + x
    elseif length(adjg) == 2
        return Polynomial([BigInt(1), BigInt(2)])
    else
        root = find_inner_vertex(adjg)
        IP_VISIT(adjg, root)
        return I[root]
    end
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::OBJ_TYPE
    return obj  # Return the original object unchanged
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

function empty_starting_point()::OBJ_TYPE
    #TODO: recover sample generation   
    prufer = rand(1:(N+2), N)
    return prufer
end
