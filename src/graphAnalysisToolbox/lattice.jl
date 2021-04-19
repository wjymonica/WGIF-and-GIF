using SparseArrays
using LinearAlgebra
function meshgrid(x, y)
    (len_x, ) = size(x)
    (len_y, ) = size(y)
    out_x = x' .* ones(Int64, len_y)
    out_y = ones(Int64, len_x)' .* y
    return out_x, out_y
end

function lattice(X, Y, connect)
    if X*Y==1
        points=[1;1]
        edges=[]
        return
    end

    if connect < 2
        rangex = 0:(X-1)
        rangey = 0:(Y-1)
        ##[x y]=meshgrid(rangey,rangex); in matlab
        x, y = meshgrid(rangey,rangex)
        points = [x[:] y[:]]
        N = X * Y

        #Connect points
        edges = [Array(1:N) (Array(1:N) .+1)]
        edges = [(edges[:,1]) edges[:,2]; Array(1:N) (Array(1:N) .+X)]
        if connect == 1 #8-connected
            border = 1:N
            border1 = findall(x->x!=0, (mod.(border, X) .-1))
            border2 = findall(x->x!=0, (mod.(border, X)))
            edges=[edges; border1 (border1 .+X .-1); border2 (border2 .+X .+1)];
        end
        excluded = findall(x->x!=0, (edges[:,1] .>N) .|(edges[:,1] .<1) .|(edges[:,2] .>N) .|(edges[:,2] .<1))
        list_excluded = [excluded; Array(X:X:((Y-1)*X))]
        edges = edges[setdiff(1:end, list_excluded), :]
    else
        #Generate points
        X = X + 2 * connect
        Y = Y + 2 * connect
        rangex = 0:(X-1)
        rangey = 0:(Y-1)
        x, y=meshgrid(rangey,rangex);
        points = [x[:] y[:]]
        N = X * Y

        #Generate adjacency matrix
        distances = (points[:,1] .- floor(Y/2)).^2 .+ (points[:,2] .- floor(X/2)).^2
        inside = findall(x->x<=(connect^2), distances)
        adjVec = sparse(ones(size(inside)), inside, 1, 1, N)
        midPoint=findall(x->x!=0, ((points[:,1] .==floor(Y/2)) .& (points[:,2] .==floor(X/2))))[1]
        #Adjust vector to be centered around diagonal
        adjVec=adjVec[:, [midPoint:N; 1:(midPoint-1)]]
        adjVec[:,1] .=0; #Remove self-connection
        adjVec = dropzeros(adjVec)
        W=circulant(adjVec[:]);

        #Remove phantom points
        index = 1:N
        tmp = 1:X:N
        tmp = tmp * ones(Int64, 1, connect) + ones(Int64, length(tmp)) * (0:(connect-1))'
        tmp2 = tmp .+ X .- connect
        list_excluded = [tmp[:]; tmp2[:]; Array(1:(X*connect)); Array(N .+(-(X*connect-1):0))]
        index = index[setdiff(1:end, list_excluded)]
        edges=adjtoedges(W[index, index])
        points=points[index, :]
    end
    return points, edges
end

function ind2subv(shape, indices)
    indx = CartesianIndices(shape)[indices]
    return hcat(getindex.(indx, 1), getindex.(indx,2))
end
#ind2subv(shape, indices) = Array(CartesianIndices(shape)[indices])

function circulant(v)
    N=length(v)
    #Build index
    index = findall(x->x!=0, v)
    diagFact = N+1
    tmp = (diagFact) .* (0:(N-1))
    masterIndex = index * ones(Int64,1,N) + ones(Int64,length(index)) * tmp';
    masterIndex[findall(x->x>(N^2), masterIndex)] .= masterIndex[1,1];
    #Build A
    i = ind2subv((N, N), masterIndex[:])[:,1]
    j = ind2subv((N, N), masterIndex[:])[:,2]
    lowerIndex = findall(x->x!=0, j.<=i)
    W = sparse(i[lowerIndex], j[lowerIndex], 1, N, N)
    W = (W+W').!=0
    #If vec is not sparse, make full again
    if ~issparse(v)
        W = Matrix(W)
    end
    return W
end

function adjtoedges(W)
    indx=findall(x->x!=0, triu(W));
    edges = hcat(getindex.(indx, 1), getindex.(indx,2))
    return edges
end


function adaptativelattice(X, Y, W)
    #W is (XY) × (XY)
    max_nei = 15
    map = sum(W, dims = 2)
    if !isapprox(norm(W),0)
        map = normalize(map, 1000)
    end
    map = reshape(map, X, Y)
    num_nei = Int.(round.(map .*max_nei))
    nei = Array{Any}(undef,max_nei+1)
    for i = 0:max_nei
        nei[i+1] = gen_nei(i)
    end
    edges = Array{Int64}(undef, 0, 2)
    wid = Y
    hei = X
    for j in 1:wid
        print("at colum", j, "\n")
        for i in 1:hei
            #print("at point ", i, ", ", j, "\n")
            edges = addedges(edges, nei[num_nei[i, j]+1], i, j, wid, hei)
        end
    end
    return edges
end

function addedges(edges, nei, i, j, wid, hei)
    nx, _ = size(nei)
    for m  = 1:nx
        kh = nei[m, 1]
        kw = nei[m, 2]
        neih = i + kh
        neiw = j + kw
        if ((neih<=hei) && (neih>=1) && (neiw>=1) && (neiw<=wid))
            edges = [edges; [(j-1)*hei+i (neiw-1)*hei+neih]]
        end
    end
    return edges
end

function gen_nei(nei)
    n = Array{Int64}(undef, 0, 2)
    if nei == 0
        n = [1 0; 0 1]
    elseif nei == 1
        n = [1 0; 0 1; 1 1; -1 1]
    elseif nei > 1
        for i= 0:nei
            for j= 0:nei
                if ((i)^2+(j)^2<=nei^2)
                    n = [n; i j]
                end
            end
        end
        for i= -nei:-1
            for j= 1:nei
                if ((i)^2+(j)^2<=nei^2)
                    n = [n; i j]
                end
            end
        end
        n = n[2:end,:]
    end
    return n
end
