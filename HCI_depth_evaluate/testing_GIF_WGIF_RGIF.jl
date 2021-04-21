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
include("./src/wgif.jl")
include("./src/rgif.jl")
include("./src/rgif_rev.jl")
include("./src/gif.jl")
include("./src/joint_BF.jl")
include("./HCI_depth_evaluate/fench_data.jl")

# ==================================
# Measurements
# sum([ abs(img - Ground_Truth) ])
# percentage of bad matching pixels
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
# initialization
BF_trr = zeros(size(tof_radial_rectified));
GIF_trr = zeros(size(tof_radial_rectified));
WGIF_trr = zeros(size(tof_radial_rectified));
RGIF_trr = zeros(size(tof_radial_rectified));

# ==================================
# GIF testing operation
# Parameter setting
Operation_name = "GIF_self_guided"
result_path = "./HCI_depth_evaluate/testing_results/GIF_self_guided_results/"
low_clim = 1
high_clim = 2.3

for r in [1,2,3,4,5,10, 15, 20]
    for e in [1/10, 1/100, 1/1000]
        # Operation
        GIF_self_guided = gif( tof_radial_rectified, tof_radial_rectified, r, e);

        # Quantitatively evaluation
        sigma = 0.1
        pixel_num, pbp = PBP(GIF_self_guided, tof_gt_radial, sigma)
        rmse = RMSE(GIF_self_guided, tof_gt_radial)

        pbp = round(pbp, digits = 4)
        rmse = round(rmse, digits = 4)

        print("=====================", "\n")
        print("Parameter setting for GIF_self_guided: r="*string(r)*", ϵ="*string(e), "\n")
        print("Bad Pixel number is ", pixel_num, ", PBP = ", pbp*100, "% \n")
        print("RMSE is ", rmse, "\n")
        print("=====================", "\n")

        # Plot
        l = @layout [a  b]
        p2 = Array{Any}(nothing, (2,1))

        p2[1] = heatmap(GIF_self_guided, title= Operation_name*"  r="*string(r)*", ϵ="*string(e), clim=(low_clim,high_clim))
        p2[2] = heatmap(GIF_self_guided - tof_gt_radial, title="Difference of GIF_self_guided with GT")

        plot(p2..., layout = l, size = (1600, 600))

        if !ispath(result_path)
            mkpath(result_path)
        end
        savefig( result_path * Operation_name * "  r="*string(r)*", ϵ="*string(e)*".png")
    end
end

# ==================================
# GIF guided by intensity testing operation
# Parameter setting
Operation_name = "GIF_intensity_guided"
result_path = "./HCI_depth_evaluate/testing_results/GIF_intensity_guided_results/"
low_clim = 1
high_clim = 2.3

for r in [1,2,3,4,5,10, 15, 20]
    for e in [1/10, 1/100, 1/1000]

        # Operation
        GIF_intens_guided = gif( tof_intensity_rectified, tof_radial_rectified, r, e);

        # Quantitatively evaluation
        sigma = 0.1
        pixel_num, pbp = PBP(GIF_intens_guided, tof_gt_radial, sigma)
        rmse = RMSE(GIF_intens_guided, tof_gt_radial)

        pbp = round(pbp, digits = 4)
        rmse = round(rmse, digits = 4)

        print("=====================", "\n")
        print("Parameter setting for GIF_intensity_guided: r="*string(r)*", ϵ="*string(e), "\n")
        print("Bad Pixel number is ", pixel_num, ", PBP = ", pbp*100, "% \n")
        print("RMSE is ", rmse, "\n")
        print("=====================", "\n")

        # Plot
        l = @layout [a  b]
        p2 = Array{Any}(nothing, (2,1))

        p2[1] = heatmap(GIF_intens_guided, title= Operation_name*"  r="*string(r)*", ϵ="*string(e), clim=(low_clim,high_clim))
        p2[2] = heatmap(GIF_intens_guided - tof_gt_radial, title="Difference of GIF_intensity_guided with GT")

        plot(p2..., layout = l, size = (1600, 600))
   
        if !ispath(result_path)
            mkpath(result_path)
        end
        savefig( result_path * Operation_name*"  r="*string(r)*", ϵ="*string(e)*".png")
    end
end

# ==================================
# Bilateral Filter testing operation
# Parameter setting
result_path = "./HCI_depth_evaluate/testing_results/BF_results/"
Operation_name = "Bilateral_Filter"
low_clim = 1
high_clim = 2.3

r_bf = 20
σ_r = 3
σ_s = 3

# Operation
# ================================
# BF for tof_radial_rectified
BF =  joint_bilateral_filter(σ_r , σ_s, r_bf, tof_radial_rectified, tof_radial_rectified);
JBF = joint_bilateral_filter(σ_r , σ_s, r_bf, tof_radial_rectified, tof_intensity_rectified);

BF = BF

# Quantitatively evaluation
sigma = 0.1
pixel_num, pbp = PBP(BF, tof_gt_radial, sigma)
rmse = RMSE(BF, tof_gt_radial)

pbp = round(pbp, digits = 4)
rmse = round(rmse, digits = 4)

print("=====================", "\n")
print("Parameter setting for Bilateral filter: r_bf="*string(r_bf)*", σ_r="*string(σ_r)*", σ_s="*string(σ_s), "\n")
print("Bad Pixel number is ", pixel_num, ", PBP = ", pbp*100, "% \n")
print("RMSE is ", rmse, "\n")
print("=====================", "\n")

# Plot
l = @layout [a  b]
p2 = Array{Any}(nothing, (2,1))

p2[1] = heatmap(BF, title= Operation_name*"  r_bf="*string(r_bf)*", σ_r="*string(σ_r)*", σ_s="*string(σ_s), 
    clim=(low_clim,high_clim))
p2[2] = heatmap(BF - tof_gt_radial, title="Difference of BF with GT", clim=(0,0.8))

plot(p2..., layout = l, size = (1600, 600))

if !ispath(result_path)
    mkpath(result_path)
end
savefig( result_path * Operation_name*"  r_bf="*string(r_bf)*", σ_r="*string(σ_r)*", σ_s="*string(σ_s)*".png")


# ==================================
# Rolling GIF testing operation
# Rolling_GIF_guided
# Parameter setting
Operation_name = "Rolling_GIF_guided"
result_path = "./HCI_depth_evaluate/testing_results/RGIF_results/"
low_clim = 1
high_clim = 2.3
step = 10

for r in [1,2,3,4,5,10, 15, 20]
    for e in [1/10, 1/100, 1/1000]
        
        # Operation
        I =  zeros(size(tof_radial_rectified));
        RGIF_trr = rgif(I, tof_radial_rectified, r, e, step);

        # Quantitatively evaluation
        sigma = 0.1
        pixel_num, pbp = PBP(RGIF_trr, tof_gt_radial, sigma)
        rmse = RMSE(RGIF_trr, tof_gt_radial)

        pbp = round(pbp, digits=4)
        rmse = round(rmse, digits=4)

        print("=====================", "\n")
        print("Parameter setting for WGIF_trr: r="*string(r)*", ϵ="*string(e)*", step="*string(step), "\n")
        print("Bad Pixel number is ", pixel_num, ", PBP = ", pbp*100, "% \n")
        print("RMSE is ", rmse, "\n")
        print("=====================", "\n")
        print("\n")

        # Plot
        l = @layout [a  b]
        p2 = Array{Any}(nothing, (2,1))

        p2[1] = heatmap(RGIF_trr, title= Operation_name*"  r="*string(r)*", ϵ="*string(e)*", step="*string(step), clim=(low_clim,high_clim))
        p2[2] = heatmap(RGIF_trr - tof_gt_radial, title="Difference of Rolling_GIF_guided with GT", clim=(0,0.8))

        plot(p2..., layout = l, size = (1600, 600))
        
        if !ispath(result_path)
            mkpath(result_path)
        end
        savefig( result_path * Operation_name*"  r="*string(r)*", ϵ="*string(e)*", step="*string(step)*".png")
    end
end

# ==================================
# Weighted GIF testing operation
# Parameter setting
Operation_name = "Weight_GIF_guided"
result_path = "./HCI_depth_evaluate/testing_results/WGIF_results/"

low_clim = 1
high_clim = 2.3

for r in [1,2,3,4,5,10, 15, 20]
    for e in [1/10, 1/100, 1/1000]

        # Operation
        WGIF_trr = wgif(tof_intensity_rectified, tof_radial_rectified, r, e);

        # Quantitatively evaluation
        sigma = 0.1
        pixel_num, pbp = PBP(WGIF_trr, tof_gt_radial, sigma)
        rmse = RMSE(WGIF_trr, tof_gt_radial)

        pbp = round(pbp, digits=4)
        rmse = round(rmse, digits=4)

        print("=====================", "\n")
        print("Parameter setting for WGIF_trr: r="*string(r)*", ϵ="*string(e), "\n")
        print("Bad Pixel number is ", pixel_num, ", PBP = ", pbp*100, "% \n")
        print("RMSE is ", rmse, "\n")
        print("=====================", "\n")

        # Plot
        l = @layout [a  b]
        p2 = Array{Any}(nothing, (2,1))

        p2[1] = heatmap(WGIF_trr, title= Operation_name*"  r="*string(r)*", ϵ="*string(e), clim=(low_clim,high_clim))
        p2[2] = heatmap(WGIF_trr - tof_gt_radial, title="Difference of Weight_GIF_guided with GT", clim=(0,0.8))

        plot(p2..., layout = l, size = (1600, 600))
        
        if !ispath(result_path)
            mkpath( result_path)
        end
        savefig( result_path * Operation_name*"  r="*string(r)*", ϵ="*string(e)*".png")
    end
end