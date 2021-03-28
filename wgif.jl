include("edgeaware.jl")
include("boxfilter.jl")
function wgif(I, p, r, lambda)
    Γ = edgeaware(I, 1)
    h, w = size(I)
    N = boxfilter(ones(h, w), r) #num of pixels in the window
    #step 1:
    mean_I = calc_mean(I, r, N)
    mean_p = calc_mean(p, r, N)
    corr_I = calc_mean(I .*I, r, N)
    corr_Ip = calc_mean(I .*p, r, N)

    #step 2:
    var_I = corr_I - mean_I .* mean_I
    cov_Ip = corr_Ip - mean_I .* mean_p

    #step 3:
    a = cov_Ip ./ (var_I .+ lambda ./ Γ)
    b = mean_p - a .* mean_I

    #step 4:
    mean_a = calc_mean(a, r, N)
    mean_b = calc_mean(b, r, N)

    #step 5:
    q = mean_a .* I + mean_b
    q = q
    return q
end

function wgif_1(I, p, r, lambda)
    Γ = edgeaware(I, 1)+edgeSobel(colorview(Gray, I))
    h, w = size(I)
    N = boxfilter(ones(h, w), r) #num of pixels in the window
    #step 1:
    mean_I = calc_mean(I, r, N)
    mean_p = calc_mean(p, r, N)
    corr_I = calc_mean(I .*I, r, N)
    corr_Ip = calc_mean(I .*p, r, N)

    #step 2:
    var_I = corr_I - mean_I .* mean_I
    cov_Ip = corr_Ip - mean_I .* mean_p

    #step 3:
    a = cov_Ip ./ (var_I .+ lambda ./ Γ)
    b = mean_p - a .* mean_I

    #step 4:
    mean_a = calc_mean(a, r, N)
    mean_b = calc_mean(b, r, N)

    #step 5:
    q = mean_a .* I + mean_b
    q = q
    return q
end

function wgif_2(I, p, r, lambda)
    Γ = edgeaware(I, 1)+edgeLap(colorview(Gray, I))
    h, w = size(I)
    N = boxfilter(ones(h, w), r) #num of pixels in the window
    #step 1:
    mean_I = calc_mean(I, r, N)
    mean_p = calc_mean(p, r, N)
    corr_I = calc_mean(I .*I, r, N)
    corr_Ip = calc_mean(I .*p, r, N)

    #step 2:
    var_I = corr_I - mean_I .* mean_I
    cov_Ip = corr_Ip - mean_I .* mean_p

    #step 3:
    a = cov_Ip ./ (var_I .+ lambda ./ Γ)
    b = mean_p - a .* mean_I

    #step 4:
    mean_a = calc_mean(a, r, N)
    mean_b = calc_mean(b, r, N)

    #step 5:
    q = mean_a .* I + mean_b
    q = q
    return q
end

function wgif_3(I, p, r, lambda)
    Γ = edgeaware(I, 1)+ edgeHPF(colorview(Gray, I))
    h, w = size(I)
    N = boxfilter(ones(h, w), r) #num of pixels in the window
    #step 1:
    mean_I = calc_mean(I, r, N)
    mean_p = calc_mean(p, r, N)
    corr_I = calc_mean(I .*I, r, N)
    corr_Ip = calc_mean(I .*p, r, N)

    #step 2:
    var_I = corr_I - mean_I .* mean_I
    cov_Ip = corr_Ip - mean_I .* mean_p

    #step 3:
    a = cov_Ip ./ (var_I .+ lambda ./ Γ)
    b = mean_p - a .* mean_I

    #step 4:
    mean_a = calc_mean(a, r, N)
    mean_b = calc_mean(b, r, N)

    #step 5:
    q = mean_a .* I + mean_b
    q = q
    return q
end

function wgif_4(I, p, r, lambda)
    Γ = edgeaware(I, 1)+edgeCanny(colorview(Gray, I))
    h, w = size(I)
    N = boxfilter(ones(h, w), r) #num of pixels in the window
    #step 1:
    mean_I = calc_mean(I, r, N)
    mean_p = calc_mean(p, r, N)
    corr_I = calc_mean(I .*I, r, N)
    corr_Ip = calc_mean(I .*p, r, N)

    #step 2:
    var_I = corr_I - mean_I .* mean_I
    cov_Ip = corr_Ip - mean_I .* mean_p

    #step 3:
    a = cov_Ip ./ (var_I .+ lambda ./ Γ)
    b = mean_p - a .* mean_I

    #step 4:
    mean_a = calc_mean(a, r, N)
    mean_b = calc_mean(b, r, N)

    #step 5:
    q = mean_a .* I + mean_b
    q = q
    return q
end
