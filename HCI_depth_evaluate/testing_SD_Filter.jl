using Images, FileIO
using LinearAlgebra
using Plots
using HDF5
using LinearAlgebra
using LaTeXStrings
using MIRTjim: jim
using Plots; 
using Statistics;
using LocalFilters;
include("./src/applications.jl")
include("./src/wgif.jl")
include("./src/rgif.jl")
include("./src/rgif_rev.jl")
include("./src/gif.jl")
include("./src/sdfilter.jl")
include("./src/ag_sdfilter.jl")
include("./HCI_depth_evaluate/fench_data.jl")

# ==================================
# Measurements
function PBP(img, gt, sigma)
    (m, n) = size(img)
    diff = abs.(img .- gt)
    
    varify(t) = t >= sigma ? 1 : 0 
    binary = varify.(diff)
    
    output = sum(binary)
    pbp = output / (m*n)
    return output, pbp
end

function RMSE(img, gt)
    output = sqrt(mean( (img .- gt).^2 ) )
    return output
end

# ==================================
# Load data
file = h5open("./HCI_depth_evaluate/dataset/hcibox_stereo_tof_fullCalibrated.h5","r")
tof_gt_radial, tof_intensity_rectified, tof_radial_rectified = Fench_data(file)

# ==================================
# GIF_trr initialization
BF_trr = zeros(size(tof_radial_rectified));
GIF_trr = zeros(size(tof_radial_rectified));
WGIF_trr = zeros(size(tof_radial_rectified));
RGIF_trr = zeros(size(tof_radial_rectified));

# ==================================
# SD_filter testing
# Parameter setting
Operation_name = "SD_Filter"
nei= 0;                # 0: 4-neighbor 1: 8-neighbor
lambda = 10;           # smoothness parameter
mu = 500;              # bw for static guidance
nu = 200;              # bw for dynamic guidance
step=20;
issp = true;

low_clim = 1; high_clim = 2.3;


u0 = ones(size(tof_radial_rectified));
SD_result = sdfilter(tof_intensity_rectified,u0,tof_radial_rectified, nei,lambda,mu,nu,step,true);

# Quantitatively evaluation
sigma = 0.1
pixel_num, pbp = PBP(SD_result, tof_gt_radial, sigma)
rmse = RMSE(SD_result, tof_gt_radial)

pbp = round(pbp, digits=4)
rmse = round(rmse, digits=4)

print("=====================", "\n")
print("Parameter setting for WGIF_trr: μ="*string(mu)*", ν=" * string(nu) * ", step="*string(step), "\n")
print("Bad Pixel number is ", pixel_num, ", PBP = ", pbp*100, "% \n")
print("RMSE is ", rmse, "\n")
print("=====================", "\n")

# Plot
l = @layout [a  b]
p2 = Array{Any}(nothing, (2,1))

p2[1] = heatmap(SD_result, title= Operation_name*"  μ="*string(mu)*", ν="*string(nu)*", step="*string(step), clim=(low_clim,high_clim))
p2[2] = heatmap(SD_result - tof_gt_radial, title="Difference of SD_result with GT", clim=(0,0.8))

plot(p2..., layout = l, size = (1600, 600))

result_path = "./HCI_depth_evaluate/testing_results/SD_results/"
if !ispath(result_path)
    mkpath(result_path)
end
savefig( result_path * Operation_name* "μ="*string(mu)*", ν="*string(nu)*", step="*string(step)*".png")


# ==================================
# AG_SD_filter testing
# Parameter setting
Operation_name = "AGSD_Filter"
r = 20
e = 1/1000
I = tof_radial_rectified
RGIF_rev_trr = rgif_rev(I, tof_intensity_rectified, 20, e, 1);
sd_itr = 10

nei = 0
lambda = 10;           # smoothness parameter
mu = 500;              # bw for static guidance
nu = 200;              # bw for dynamic guidance
issp = true;
low_clim = 1; high_clim = 2.3;
u0 = ones(size(tof_radial_rectified))

# ag_sdfilter(g,u0,f,nei,lambda,sigma_g,sigma_u,itr,issparse, r=20, e=1/1000, num_rolloing=1)
sd_itr = 10
AGSD_result = sdfilter(RGIF_rev_trr, u0, tof_radial_rectified, nei, lambda, mu, nu, sd_itr, issp);

# Quantitatively evaluation
sigma = 0.1
pixel_num, pbp = PBP(AGSD_result, tof_gt_radial, sigma)
rmse = RMSE(AGSD_result, tof_gt_radial)

pbp = round(pbp, digits=4)
rmse = round(rmse, digits=4)

print("=====================", "\n")
print("Parameter setting for WGIF_trr: r="*string(r)*", ϵ="*string(e)*", step="*string(step), "\n")
print("Bad Pixel number is ", pixel_num, ", PBP = ", pbp*100, "% \n")
print("RMSE is ", rmse, "\n")
print("=====================", "\n")

# Plot
l = @layout [a  b]
p2 = Array{Any}(nothing, (2,1))

p2[1] = heatmap(AGSD_result, title= Operation_name*"  r="*string(r)*", ϵ="*string(e)*", step="*string(step), clim=(low_clim,high_clim))
p2[2] = heatmap(AGSD_result - tof_gt_radial, title="Difference of AGSD_result with GT", clim=(0,0.8))

plot(p2..., layout = l, size = (1600, 600))

result_path = "./HCI_depth_evaluate/testing_results/AGSD_results/"
if !ispath(result_path)
    mkpath(result_path)
end
savefig( result_path * Operation_name*"  r="*string(r)*", ϵ="*string(e)*", step="*string(step)*".png")