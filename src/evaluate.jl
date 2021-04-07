include("gif.jl")
include("rgif.jl")
include("wgif.jl")
using Noise
"""
step_image = zeros(1024, 1024)
step_image[:, 1:512] .= 1
noisy_image = add_gauss(step_image, 0.1)

test1 = wgif_1(noisy_image, noisy_image, 10, 1/10)
test2 = wgif_2(noisy_image, noisy_image, 10, 1/10)
test3 = wgif_3(noisy_image, noisy_image, 10, 1/10)
test4 = wgif_4(noisy_image, noisy_image, 10, 1/10)
#colorview(Gray, test1)
#range = 250:750
line1 = [step_image[513, :] noisy_image[513, :] test[513, :] test1[513, :]]
line2 = [step_image[513, :] noisy_image[513, :] test[513, :] test2[513, :]]
line3 = [step_image[513, :] noisy_image[513, :] test[513, :] test3[513, :]]
line4 = [step_image[513, :] noisy_image[513, :] test[513, :] test4[513, :]]
lines = [step_image[513, :] noisy_image[513, :] test1[513, :] test2[513, :] test3[513, :] test4[513, :]]
label = ["Step function" "Noisy signal" "WGIF(Sobel)" "WGIF(Laplacian)" "WGIF(HPF)" "WGIF(Canny)"]
label1 = ["Step function" "Noisy signal" "WGIF" "WGIF(Sobel)"]
label2 = ["Step function" "Noisy signal" "WGIF" "WGIF(Laplacian)"]
label3 = ["Step function" "Noisy signal" "WGIF" "WGIF(HPF)"]
label4 = ["Step function" "Noisy signal" "WGIF" "WGIF(Canny)"]
#lines = [step_image[5, range] noisy_image[5, range] test1[5, range] test2[5, range] test3[5, range]]
#p1 = plot(lines, label = ["Step function" "Noisy signal" "GIF" "WGIF" "WGIF-m"], xlabel = "Pixel coordinate at the horizontal direction", ylabel = "Pixel intensity")
#p1 = plot(lines, label = ["Step function" "Noisy signal" "GIF"], xlabel = "Pixel coordinate at the horizontal direction", ylabel = "Pixel intensity")

p = Array{Any}(undef, 3, 2)
p[1,1] = jim(transpose(step_image), "Original image")
p[2,1] = jim(transpose(noisy_image), "Noisy image")
p[3,1] = jim(transpose(test1), "WGIF(Sobel)")
p[1,2] = jim(transpose(test2), "WGIF(Laplacian)")
p[2,2] = jim(transpose(test3), "WGIF(HPF)")
p[3,2] = jim(transpose(test4), "WGIF(Canny)")

plot(p..., layout = (3, 2), colorbar=false)
#@time begin
#  #code
#end
#savefig("rgif4.png")
"""
I3 = load("WGIF-Julia/img/house/A.jpg")
p3 = load("WGIF-Julia/img/house/B.jpg")
I1 = Gray.(I3)
p1 = Gray.(p3)
I3 = Float64.(channelview(I3))
I1 = Float64.(I1)
p3 = Float64.(channelview(p3))
p1 = Float64.(p1)

@time begin
  wgif_4(I1, p1, 1, 1/100)
end

@time begin
  wgif_4(I3[1,:,:], p3[1,:,:], 1, 1/100)
  wgif_4(I3[2,:,:], p3[2,:,:], 1, 1/100)
  wgif_4(I3[3,:,:], p3[3,:,:], 1, 1/100)
end
@time begin
  wgif_4(I1, I1, 1, 1/100)
end
@time begin
  wgif_4(I3[1,:,:], I3[1,:,:], 1, 1/100)
  wgif_4(I3[2,:,:], I3[2,:,:], 1, 1/100)
  wgif_4(I3[3,:,:], I3[3,:,:], 1, 1/100)
end
