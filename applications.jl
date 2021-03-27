include("wgif.jl")
include("gif.jl")
using Images
using PaddedViews

function detail_enhancement(I, p, r, lambda_e, θ, method)
    c, h, w = size(I)
    G = zeros(c, h, w)
    for i = 1:c
        if method == "wgif"
            G[i, :, :] = wgif(I[i, :, :], p[i, :, :], r, lambda_e)
        elseif method == "gif"
            G[i, :, :] = gif(I[i, :, :], p[i, :, :], r, lambda_e)
        else
            G[i, :, :] = G[:, :, i]
        end
    end
    out_I = (I-G) * θ + G
end

function getAirLight(I_color, I_gray, ratio)
    c, h, w,  = size(I_color)
    I_b = I_gray .> (ratio * maximum(I_gray))
    count = sum(I_b)
    out = zeros(c)
    for i = 1:c
        sum_b = sum(I_color[i,:,:] .* I_b)
        out[i] = sum_b / count
    end
    return out
end

function MinFilter(I, r)
    h, w = size(I);
    out_I = zeros(h, w);
    for i = 1:h
        for j = 1:w
            x1 = max(min(i-r, h), 1)
            x2 = max(min(i+r, h), 1)
            y1 = max(min(j-r, w), 1)
            y2 = max(min(j+r, w), 1)
            out_I[i, j] = minimum(I[x1:x2, y1:y2]);
        end
    end
    return out_I
end

function dehazing(I, method, r = 10, e = 7, ratio = 0.98, omega = 31/32, lambda = 1/1000, haze_level = 0.03125)
    c, h, w = size(I)
    I_rgb = Float64.(I) ./ 255;
    I_min = mapslices(minimum, I_rgb; dims = [1])[1, :, :]
    I_d = MinFilter(I_min, e)
    ac = getAirLight(I_rgb, I_d, 0.98)
    t = 1 .- omega * (I_d ./ minimum(ac))
    I_gray = Float64.(Gray.(colorview(RGB,I_rgb)))
    t_smoothed = t
    if method == "wgif"
        t_smoothed = wgif(t, I_gray, r, lambda)
    end
    if method == "gif"
        t_smoothed = wgif(t, I_gray, r, lambda)
    end
    t_adjusted = t_smoothed .^(1+haze_level)
    amplify_factor0 = max.(0.1, 1 ./ t_adjusted .- 1)
    amplify_factor = zeros(3, size(amplify_factor0)[1],size(amplify_factor0)[2])
    amplify_factor[1, :, :] = amplify_factor0
    amplify_factor[2, :, :] = amplify_factor0
    amplify_factor[3, :, :] = amplify_factor0
    ac_rep = zeros(3, h, w)
    ac_rep[1, :, :] .= ac[1]
    ac_rep[2, :, :] .= ac[2]
    ac_rep[3, :, :] .= ac[3]
    G = I_rgb + amplify_factor .* (I_rgb - ac_rep)
    return G.*255
end

#function exposure_fusion()
