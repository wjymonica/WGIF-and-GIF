using LinearAlgebra
function makeweightsL2(edges,vals,valScale,points=zeros(size(edges,1),1),geomScale=0,EPSILON=1e-5)
    if valScale > 0
        valDistance = sum((vals[edges[:,1],:] .- vals[edges[:,2],:]).^2, dims=2)
        if !isapprox(norm(valDistance), 0)
            valDistance = normalize(valDistance, 1000)
        end
    else
        valDistance = zeros(size(edges,1), 1)
        valScale = 0
    end

    geomDistances=zeros(size(edges,1),1);
    geomScale=0;

    weights=exp.(-(geomScale * geomDistances .+ valScale * valDistance)) .+ EPSILON

    return weights
end
