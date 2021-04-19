using LinearAlgebra
using SparseArrays
function adjacency(edges, weights, N)
    #if N < 0
    #    N=maximum(edges)
    #end

    #W = Array{Float64, 2}(undef, N, N)
    #W[:] .= 0
    #for i = 1:N
    #    W[edges[i,1], edges[i,2]] = weights[i]
    #    W[edges[i,2], edges[i,1]] = weights[i]
    #end
    W = sparse([edges[:, 1]; edges[:, 2]], [edges[:, 2]; edges[:, 1]], [weights; weights][:], N, N)

    return W

end
