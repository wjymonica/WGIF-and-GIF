include("gif.jl")
# when self-guided: use rg = rgif(img, zeros(size(img)),3,0.5,4)
function rgif(I, p, r, e,iteration=4)           
    res = I;
    
    for i in 1:iteration
        print("Rolling GIF iteration ", i, "...");         
        res = gif(res, p, r, e);       
    end    

    return res
end