include("gif.jl")
# difference from rgif in that it make guidance fixed and keep
# ``input'' image updating
function rgif_rev(I, p, r, e,iteration=4)           
    res = p;
    
    for i in 1:iteration
        print("Reverse Rolling GIF iteration ", i, "...");         
        res = gif(I, res, r, e);       
    end    

    return res
end