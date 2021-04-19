include("gif.jl")
include("rgif.jl")
include("wgif.jl")
include("sdfilter.jl")
include("ag_sdfilter.jl")
using Noise
using MIRTjim:jim
using Plots
using Images, FileIO
using Statistics
using TestImages
using ImageBinarization
using LinearAlgebra
using ImageMagick

nei=5;                # 0: 4-neighbor 1: 8-neighbor
lambda = 10;           # smoothness parameter
mu = 500;              # bw for static guidance
nu = 200;              # bw for dynamic guidance
step=20;
issp = true;
#m, n, = size(img)


img =  ImageMagick.load("/Users/wangjingying/Desktop/EECS556/sdfilter/imgs/depthSR/books.bmp")

img = Float64.(Gray.(img))
noisy_image = add_gauss(img, 0.1)
u0 = ones(size(noisy_image))
colorview(Gray, img)

out_gif = gif(noisy_image, noisy_image, 10, 1/10)
out_wgif = wgif(noisy_image, noisy_image, 10, 1/10)
out_wgif1 = wgif_1(noisy_image, noisy_image, 10, 1/10)
out_wgif2 = wgif_2(noisy_image, noisy_image, 10, 1/10)
out_wgif3 = wgif_3(noisy_image, noisy_image, 10, 1/10)
out_wgif4 = wgif_4(noisy_image, noisy_image, 10, 1/10)
out_sd = sdfilter(noisy_image,u0,noisy_image,nei,lambda,mu,nu,20,issp)
out_asd = asd(noisy_image,u0,noisy_image,nei,lambda,mu,nu,20,issp)
out_agsd = agsd(noisy_image,u0,noisy_image,nei,lambda,mu,nu,20,issp)
