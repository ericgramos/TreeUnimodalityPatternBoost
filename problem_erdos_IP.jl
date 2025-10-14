# This file takes a prufer code and prints the independence polynomial, reward, where failure happens, 
# and log-concavity metric at all places in the independence sequence. All computations are done with BigInt, just to be safe.

include("constants.jl")
using JSON
using Polynomials
using DataStructures
using Random

#Code is checking for breakage at index ratio*n - offset
const ratio = 1/2
const offset = 0

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


function IP(adjg::Vector{Vector{Int}})::Tuple{Polynomial{BigInt}, Integer, Integer}
    x = Polynomial{Int}([BigInt(0), BigInt(1)]) 
    I = Dict{Int, Polynomial{BigInt}}()
    V = Dict{Int, Polynomial{BigInt}}()
    N = Dict{Int, Polynomial{BigInt}}()
    C = [BigInt[] for _ in 1:length(adjg)]
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
            V[root] = Polynomial{BigInt}([BigInt(1)])  
            N[root] = Polynomial{BigInt}([BigInt(1)])
            I[root] = Polynomial{BigInt}([BigInt(1)]) + x
            return
        end
        for u in adjg[root]
            if u ∉ visited
                push!(C[root], u)
                IP_VISIT(adjg, u)
            end
        end

        left = Polynomial{BigInt}([BigInt(1)])
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
        return Polynomial{BigInt}([BigInt(1)]) + x
    elseif length(adjg) == 2
        return Polynomial{BigInt}([BigInt(1), BigInt(2)])
    else
        root = find_inner_vertex(adjg)
        IP_VISIT(adjg, root)
        alpha = degree(I[root])
        coeffs_array = coeffs(I[root])
        index = 1
        while coeffs_array[index] <= coeffs_array[index+1]
            index += 1
        end
        mode_index = index
        # mode_index is i+1 because the sequence starts at 0
        #mode_index = findmax(coeffs_array)[2]
        return I[root], alpha, mode_index
    end
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    prufer = obj
    adjg = prufer_to_adjacency_list(prufer)
    n = length(adjg)

    poly, alpha, mode_index = IP(adjg)
    coeffs_array = BigInt.(coeffs(poly))
    #TODO
    # checking at n/2 floor index (sequence starts at 0)
    index = floor(Int, ratio*(N+2))-offset+1 
    
    # Pad coeffs_array with zeros if needed
    expected_length = index + 1
    if length(coeffs_array) < expected_length
        append!(coeffs_array, zeros(Int, expected_length - length(coeffs_array)))
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

#Enter the prufer code of the tree you are testing here:
prufer = [18,21,4,22,7,20,16,10,20,24,17,6,25,25,20,25,20,19,1,22,25,19,22,13]

N = length(prufer) 
adjg = prufer_to_adjacency_list(prufer)
poly, alpha, mode_index = IP(adjg)
coeffs_array = BigInt.(coeffs(poly))
reward = reward_calc(prufer)
log_count = []
for i in 2: alpha
    index = i+1 
    expected_length = index + 1
    if length(coeffs_array) < expected_length
        append!(coeffs_array, zeros(Int, expected_length - length(coeffs_array)))
    end
    a1 = coeffs_array[index-1]
    a2 = coeffs_array[index]
    a3 = coeffs_array[index+1]
    push!(log_count, (a1*a3)-(a2^2))
end

println("prufer: ", prufer)
println("poly: ", poly)
println(log_count)
println("checking at: ", floor(Int, ratio*(N+2)-offset))
println("reward: ", reward)