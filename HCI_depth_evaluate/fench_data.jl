function Fench_data(file)
    # ================================
    # Open file.h5
    # println("names n",names(file))

    # ================================
    # Read datasets
    tof_gt_depth = read(file, "tof_groundTruthDepth");
    tof_gt_radial = read(file, "tof_groundTruthRadialDistance");
    tof_intensity_rectified = read(file, "tof_intensity_rectified"); # 200×200×1×1×1 
    tof_radial_rectified = read(file, "tof_radialDistance_rectified"); # 200×200×1×1×1 

    tof_gt_depth = reshape(tof_gt_depth, (200, 200));
    tof_gt_radial = reshape(tof_gt_radial, (200, 200));
    tof_intensity_rectified = reshape(tof_intensity_rectified, (200, 200));
    tof_radial_rectified = reshape(tof_radial_rectified, (200, 200));

    # ================================
    # Crop the matrix margin
    tof_gt_depth = tof_gt_depth[16:(16+162), 18:(18+162)];  # depth GT
    tof_gt_radial = tof_gt_radial[16:(16+162), 18:(18+162)];  # radius depth GT
    tof_intensity_rectified = tof_intensity_rectified[16:(16+162), 18:(18+162)];  # intensity
    tof_radial_rectified = tof_radial_rectified[16:(16+162), 18:(18+162)];  # depth

    # ================================
    # Modify to shift the bias of tof_radial_rectified
    mean_gt_radial = mean(tof_gt_radial);
    mean_radial_rectified = mean(tof_radial_rectified);
    tof_radial_rectified = tof_radial_rectified .- mean_radial_rectified .+ mean_gt_radial;

    # ==================================
    # rescale the intensity image
    max_intensity = maximum(tof_intensity_rectified)
    min_intensity = minimum(tof_intensity_rectified)

    max_depth = maximum(tof_radial_rectified)
    min_depth = minimum(tof_radial_rectified)

    alpha = (max_depth - min_depth) / (max_intensity - min_intensity)
    beta = min_depth - alpha * min_intensity
    tof_intensity_rectified = alpha .* tof_intensity_rectified .+ beta

    return tof_gt_radial, tof_intensity_rectified, tof_radial_rectified
end