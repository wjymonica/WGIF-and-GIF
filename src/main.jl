using Images, FileIO
using Statistics
using TestImages
using ImageBinarization
using LinearAlgebra
using Plots
include("applications.jl")
function rgb2mat_256(path)
    img = load(path)
    img = channelview(img) .*255
    return img
end

r = 60
e = 7
ratio = 0.98
omega = 31/32
lambda = 1/1000
haze_level = 0.03125

img = load("WGIF-Julia/img/tulips.bmp")
img = Float64.(channelview(img))

#img[1, :, :] = wgif_m(img[1, :, :], img[1, :, :], 10, 1/1000)
#img[2, :, :] = wgif_m(img[2, :, :], img[2, :, :], 10, 1/1000)
#img[3, :, :] = wgif_m(img[3, :, :], img[3, :, :], 10, 1/1000)
G1 = detail_enhancement(img, img, 16, 1/128, 4, "wgif")
#G2 = dehazing(img, "gif")
colorview(RGB, G1)
