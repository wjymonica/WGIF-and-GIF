include("colorspace.jl")
include("graphAnalysisToolBox/adjacency.jl")
include("graphAnalysisToolBox/lattice.jl")
include("graphAnalysisToolBox/makeweightsL2.jl")

function Sparse(vec1, vec2, A)
    W = Array{typeof(vec1[1]), 2}(undef, size(vec1)[1], size(vec2)[1])
    W = A
    W = sparse(W)
    return W
end


function sdfilter(g,u0,f,nei,lambda,sigma_g,sigma_u,itr,issparse)

    _, _ = size(g);
    X, Y = size(u0);

    N = X .* Y;
    _, edges = lattice(X,Y,nei);

    fVals = reshape(f, N);
    if (issparse)
        #A = zeros(N)
        #fval_ind = findall(x->x>0, fVals)
        #A[fval_ind] = 1;
        C=sparse(1:N, 1:N, ones(N));
        F = C * float64.(fVals);
    else
        C = Sparse(1:N,1:N,ones(N));
        F = float64.(fVals);
    end

    gVals = reshape(g, N);

    weights_g = makeweightsL2(edges,gVals,sigma_g);

    gW = adjacency(edges,weights_g,N);

    #edges = adaptativelattice(X, Y, gW)

    #gVals = reshape(g, N);

    #weights_g = makeweightsL2(edges,gVals,sigma_g);

    #gW = adjacency(edges,weights_g,N);


    print(1, "lambda: ", lambda, ", # of steps: ", itr, "\n")
    for i=1:itr
        print(1, i, itr, "\n")
        #pause(.1)


        uVals = reshape(u0, N);
        weights_u = makeweightsL2(edges,uVals,sigma_u);
        uW = adjacency(edges,weights_u,N);

        W = gW.*uW;
        D  = sparse(1:N, 1:N, sum(W, dims = 2)[:]);
        L = D-W;

        R = (C .+ lambda.*L);
        U = (R \ F);

        u = reshape(U, X, Y);
        u0 = u;
    end

    print('\n');
    return u0
end

function asd(g,u0,f,nei,lambda,sigma_g,sigma_u,itr,issparse)
    _, _ = size(g);
    X, Y = size(u0);

    N = X .* Y;
    _, edges = lattice(X,Y,nei);

    fVals = reshape(f, N);
    if (issparse)
        #A = zeros(N)
        #fval_ind = findall(x->x>0, fVals)
        #A[fval_ind] = 1;
        C=sparse(1:N, 1:N, ones(N));
        F = C * float64.(fVals);
    else
        C = Sparse(1:N,1:N,ones(N));
        F = float64.(fVals);
    end
    print("Original edge size: ", size(edges), "\n")
    gVals = reshape(g, N);

    weights_g = makeweightsL2(edges,gVals,sigma_g);

    gW = adjacency(edges,weights_g,N);

    edges = adaptativelattice(X, Y, gW)
    print("Adaptative edge size: ", size(edges), "\n")

    gVals = reshape(g, N);

    weights_g = makeweightsL2(edges,gVals,sigma_g);

    gW = adjacency(edges,weights_g,N);


    print(1, "lambda: ", lambda, ", # of steps: ", itr, "\n")
    for i=1:itr
        print(1, i, itr, "\n")
        #pause(.1)


        uVals = reshape(u0, N);
        weights_u = makeweightsL2(edges,uVals,sigma_u);
        uW = adjacency(edges,weights_u,N);

        W = gW.*uW;
        D  = sparse(1:N, 1:N, sum(W, dims = 2)[:]);
        L = D-W;

        R = (C .+ lambda.*L);
        U = (R \ F);

        u = reshape(U, X, Y);
        u0 = u;
    end

    print('\n');
    return u0
end

function agsd(g,u0,f,nei,lambda,sigma_g,sigma_u,itr,issparse, r=20, e=1/1000, num_rolloing=1)
    # g: guidance image
    # u0: initial estimat
    # nei: window size
    # lambda: a regularization parameter
    # sigma_g: bandwidth parameter for g
    # sigma_u: a bandwidth parameter for u
    # itr: number of iterations
    # issparse: true for downsampling; false for usual case
    # r: radius (same in GIF)
    # e: regularization (sam in GIF)
    # num_rolloing: number of iteration of reverse RGIF


    g = rgif_rev(f, g, r, e, num_rolloing)

    _, _ = size(g);
    X, Y = size(u0);

    N = X .* Y;
    _, edges = lattice(X,Y,nei);

    fVals = reshape(f, N);
    if (issparse)
        #A = zeros(N)
        #fval_ind = findall(x->x>0, fVals)
        #A[fval_ind] = 1;
        C=sparse(1:N, 1:N, ones(N));
        F = C * float64.(fVals);
    else
        C = Sparse(1:N,1:N,ones(N));
        F = float64.(fVals);
    end

    gVals = reshape(g, N);

    weights_g = makeweightsL2(edges,gVals,sigma_g);

    gW = adjacency(edges,weights_g,N);


#    print("lambda: ", lambda, ", # of steps: ", itr, "\n")
    for i=1:itr
#        print("iteration # ", i, "\n")
        #pause(.1)

        uVals = reshape(u0, N);
        weights_u = makeweightsL2(edges,uVals,sigma_u);
        uW = adjacency(edges,weights_u,N);

        W = gW.*uW;
        D  = sparse(1:N, 1:N, sum(W, dims = 2)[:]);
        L = D-W;

        R = (C .+ lambda.*L);
        U = (R \ F);

        u = reshape(U, X, Y);
        u0 = u;
    end

    print('\n');
    return u0
end




# # --------------------------------------
# # This is for application about conversion of Image to a MAtrix{Float64}
# nei= 0;
# lambda = 10;
# sigma_g = 500;
# nu = 200;
# itr=20;
# issparse = true;


# img_name = "./imgs/depthSR/books.bmp"
# dep_name = "./imgs/depthSR/books_depth.bmp"

# img = float64.(ImageMagick.load(img_name))
# dep = float64.(ImageMagick.load(dep_name))

# img = Float64.(Gray.(img[:,:]))
# dep = Float64.(Gray.(dep[:,:]))

# m, n = size(img)
# u0=ones(m,n);
# #-----------------------------------------
