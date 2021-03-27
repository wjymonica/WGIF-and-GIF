using Images, FileIO
using LinearAlgebra
using Plots
include("applications.jl")
function rgb2mat_256(path)
    img = load(path)
    img = channelview(img) .*255
    return img
end

r = 10
e = 7
ratio = 0.98
omega = 31/32
lambda = 1/1000
haze_level = 0.03125

img = load("./img/forest.jpg")
img = channelview(img) .*255

G1 = detail_enhancement(img, img, 12, 1/28, 4, "wgif")
G2 = dehazing(img, "wgif")
colorview(RGB, G2./255)
