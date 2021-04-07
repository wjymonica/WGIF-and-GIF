using Images, FileIO
using Statistics
using TestImages
using ImageBinarization
using DSP
include("gif.jl")
include("boxfilter.jl")

function edgeaware(I, r)
    h, w = size(I)
    N = boxfilter(ones(h, w), r) #num of pixels in the window
    L = maximum(I)-minimum(I)
    e = (0.001*L)^2
    mean_I = calc_mean(I, r, N)
    corr_I = calc_mean(I .* I, r, N)
    var_I = corr_I - mean_I .* mean_I
    Γ = (var_I .+ e) .* sum(1 ./ (var_I .+ e)) ./ (h*w)
    return Γ
end

function edgeSobel(img)
    imgs = 64 * imfilter(Gray.(img), Kernel.sobel())
    x1 = minimum(Float64.(imgs))
    x2 = maximum(Float64.(imgs))
    return (Float64.(imgs) .- x1) ./ (x2-x1)
end

function edgeLap(img)
    imgl = imfilter(img, Kernel.Laplacian())
    imgl = Gray.(imgl)
    x1 = minimum(Float64.(imgl))
    x2 = maximum(Float64.(imgl))
    return (Float64.(imgl) .- x1) ./ (x2-x1)
end

function edgeHPF(img)
    imgh = Float64.(Gray.(img))
    h=[-1 1]
    imgh = conv(imgh,h)[:, 2:end]
    x1 = minimum(imgh)
    x2 = maximum(imgh)
    return (imgh .- x1) ./ (x2 - x1)
end

function edgeCanny(img)
    imgc = canny(Gray.(img), (Percentile(80), Percentile(20)))
    return imgc
end
