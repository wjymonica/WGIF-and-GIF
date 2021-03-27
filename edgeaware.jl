include("gif.jl")
include("boxfilter.jl")
function edgeaware(I, r)
    h, w = size(I)
    N = boxfilter(ones(h, w), r) #num of pixels in the window
    L = maximum(I)-minimum(I)
    e = (0.001*L)^2;
    mean_I = calc_mean(I, r, N)
    corr_I = calc_mean(I .*I, r, N)
    var_I = corr_I - mean_I .* mean_I
    Γ = (var_I .+ e) .* sum(1 ./ (var_I .+ e)) ./ (h*w)
    return Γ
end
