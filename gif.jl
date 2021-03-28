using LinearAlgebra
include("boxfilter.jl")
function calc_mean(I, r, N)
    return boxfilter(I, r) ./ N
end

function gif(I, p, r, e)
    I = Float64.(I)
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
    a = cov_Ip ./ (var_I .+ e)
    b = mean_p - a .* mean_I

    #step 4:
    mean_a = calc_mean(a, r, N)
    mean_b = calc_mean(b, r, N)

    #step 5:
    q = mean_a .* I + mean_b
    return q
end
