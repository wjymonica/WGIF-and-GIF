using Images, FileIO
using LinearAlgebra
using Plots
using HDF5
using LinearAlgebra
using LaTeXStrings
using MIRT: jim
using Plots; 
using Statistics;
using LocalFilters;
include("applications.jl")
include("wgif.jl")

file = h5open("./dataset/hcibox_stereo_tof_fullCalibrated.h5","r")

# println("names n",names(file))

# Read datasets
tof_radial_rectified = read(file, "tof_radialDistance_rectified"); # 200×200×1×1×1 
tof_gt_depth = read(file, "tof_groundTruthDepth");
tof_gt_radial = read(file, "tof_groundTruthRadialDistance");

tof_radial_rectified = reshape(tof_radial_rectified, (200, 200));
tof_gt_depth = reshape(tof_gt_depth, (200, 200));
tof_gt_radial = reshape(tof_gt_radial, (200, 200));

# Crop the matrix margin
tof_radial_rectified = tof_radial_rectified[16:(16+161), 18:(18+162)];
tof_gt_depth = tof_gt_depth[16:(16+161), 18:(18+162)];
tof_gt_radial = tof_gt_radial[16:(16+161), 18:(18+162)];

# Modify to shift the bias
mean_gt_radial = mean(tof_gt_radial)
mean_radial_rectified = mean(tof_radial_rectified);
tof_radial_rectified = tof_radial_rectified .- mean_radial_rectified .+ mean_gt_radial;

# BF for tof_radial_rectified
box = RectangularBox{2}(5)
BF_trr =  bilateralfilter( tof_radial_rectified, 4, 3, box)
# GIF for tof_radial_rectified
GIF_trr = wgif(tof_radial_rectified, tof_radial_rectified, 12, 1/28)

# RMSE of BF trr and tof_gt_radial
RMSE_BF = sqrt(mean(abs.(BF_trr .- tof_gt_radial).^2))
println("The RMSE of bilateral Filter is ", RMSE_BF)
# RMSE of GIF trr and tof_gt_radial
RMSE_GIF = sqrt(mean(abs.(GIF_trr .- tof_gt_radial).^2))
println("The RMSE of Guided Image Filter is ", RMSE_GIF)


pl = Array{Any}(nothing, (3,1))
pl[1] = jim(tof_radial_rectified, title="tof_radial_rectified", clim=(0,2.2))
pl[2] = jim(tof_gt_depth, title="tof_gt_depth", clim=(0,2.2))
pl[3] = jim(tof_gt_radial, title="tof_gt_radial", clim=(0,2.2))

plot(pl...)

p2 = Array{Any}(nothing, (4,1))
p2[1] = jim(tof_radial_rectified, title="tof_radial_rectified", clim=(0,2.2))
p2[2] = jim(tof_gt_radial, title="tof_gt_radial", clim=(0,2.2))
p2[3] = jim(BF_trr, title="BF trr", clim=(0,2.2))
p2[4] = jim(GIF_trr, title="GIF trr", clim=(0,2.2))

plot(p2...)