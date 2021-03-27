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
include("rgif.jl")

# ================================
# Open file.h5
file = h5open("./dataset/hcibox_stereo_tof_fullCalibrated.h5","r")
# println("names n",names(file))


# ================================
# Read datasets
tof_radial_rectified = read(file, "tof_radialDistance_rectified"); # 200×200×1×1×1 
tof_gt_depth = read(file, "tof_groundTruthDepth");
tof_gt_radial = read(file, "tof_groundTruthRadialDistance");

tof_radial_rectified = reshape(tof_radial_rectified, (200, 200));
tof_gt_depth = reshape(tof_gt_depth, (200, 200));
tof_gt_radial = reshape(tof_gt_radial, (200, 200));


# ================================
# Crop the matrix margin
tof_radial_rectified = tof_radial_rectified[16:(16+162), 18:(18+162)];
tof_gt_depth = tof_gt_depth[16:(16+162), 18:(18+162)];
tof_gt_radial = tof_gt_radial[16:(16+162), 18:(18+162)];


# ================================
# Modify to shift the bias
mean_gt_radial = mean(tof_gt_radial);
mean_radial_rectified = mean(tof_radial_rectified);
tof_radial_rectified = tof_radial_rectified .- mean_radial_rectified .+ mean_gt_radial;


# ================================
# BF for tof_radial_rectified
box = RectangularBox{2}(5);
BF_trr =  bilateralfilter( tof_radial_rectified, 4, 3, box);

# GIF for tof_radial_rectified
GIF_trr = wgif(tof_radial_rectified, tof_radial_rectified, 12, 1/28);

# Rolling GIF
r = 12;
e = 1/28;
I = ones(size(tof_radial_rectified));
RGIF_trr = rgif(I, tof_radial_rectified, r, e);


# ================================
# RMSE of BF trr and tof_gt_radial
RMSE_BF = sqrt(mean((BF_trr .- tof_gt_radial).^2))
# RMSE of GIF trr and tof_gt_radial
RMSE_GIF = sqrt(mean((GIF_trr .- tof_gt_radial).^2))
# RMSE of RGIF_trr trr and tof_gt_radial
RMSE_RGIF = sqrt(mean((RGIF_trr .- tof_gt_radial).^2))

println("The RMSE of Guided Image Filter is ", RMSE_GIF)
println("The RMSE of bilateral Filter is ", RMSE_BF)
println("The RMSE of Rolling Guided Image Filter is ", RMSE_RGIF)


# ================================
# Plot
pl = Array{Any}(nothing, (3,1))
# pl[1] = jim(tof_radial_rectified, title="tof_radial_rectified", clim=(0,2.2))
# pl[2] = jim(tof_gt_depth, title="tof_gt_depth", clim=(0,2.2))
# pl[3] = jim(tof_gt_radial, title="tof_gt_radial", clim=(0,2.2))

pl[1] = heatmap(tof_radial_rectified, title="tof_radial_rectified", clim=(1,2.2))
pl[2] = heatmap(tof_gt_depth, title="tof_gt_depth", clim=(1,2.2))
pl[3] = heatmap(tof_gt_radial, title="tof_gt_radial", clim=(1,2.2))

plot(pl...)

l = @layout [a  b  ; c  d ; t]
p2 = Array{Any}(nothing, (5,1))
# p2[1] = jim(tof_radial_rectified, title="tof_radial_rectified", clim=(0,2.2))
# p2[2] = jim(tof_gt_radial, title="tof_gt_radial", clim=(0,2.2))
# p2[3] = jim(BF_trr, title="BF trr", clim=(0,2.2))
# p2[4] = jim(GIF_trr, title="GIF trr", clim=(0,2.2))
# p2[5] = jim(RGIF_trr, title="RGIF trr", clim=(0,2.2))
p2[1] = heatmap(tof_radial_rectified, title="tof_radial_rectified", clim=(1,2.2))
p2[2] = heatmap(tof_gt_radial, title="tof_gt_radial", clim=(1,2.2))
p2[3] = heatmap(BF_trr, title="BF trr", clim=(1,2.2))
p2[4] = heatmap(GIF_trr, title="GIF trr", clim=(1,2.2))
p2[5] = heatmap(RGIF_trr, title="RGIF trr", clim=(1,2.2))

plot(p2..., layout = 5, size = (1500, 720))