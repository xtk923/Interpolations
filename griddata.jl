using LinearAlgebra
"""
griddata with linear interpolation which works similar to its matlab counterpart

`griddata(x::Array{Real, 1}, y::Array{Real, 1}, v::Array{Real, 1}, xq::Array{Real, 1}, yq::Array{Real, 1})`
"""
function griddata(x, y, v,
                  xq, yq)
    if !(length(x) == length(y) && length(y) == length(v))
        error("x, y and v must be in same dimension.")
    end
    if !(length(xq) == length(yq))
        error("xq and yq must be in same dimension.")
    end

    res = []
    for i = 1:length(xq)
        psGiven = [x y]
        dists = ((xq[i] .- x).^2 + (yq[i] .- y).^2).^0.5
        # find the nearest three by (x,y)
        idxNearest = []
        while length(idxNearest) < 3
            tmp = minimum(filter(x->!isnan(x), dists))
            append!(idxNearest, indexin(tmp, dists))
            dists[idxNearest[end]] = NaN
        end
        vertices = nothing
        for idx in idxNearest
            vertices = vertices == nothing ?
                [x[idx] y[idx] v[idx]] : vcat(vertices, [x[idx] y[idx] v[idx]])
        end
        # with three nearest points, get the expression of the plane,
        # then find the z value of f(xq[i], yq[i])
        D = det(vertices)
        d = 1
        a = -d/D * det(
            hcat([1,1,1], vertices[:, 2:3])
        )
        b = -d/D * det(
            hcat(vertices[:, 1], [1,1,1], vertices[:, 3])
        )
        c = -d/D * det(
            hcat(vertices[:, 1:2], [1,1,1])
        )
        vq = (a*xq[i] + b*yq[i] + d) / -c
        append!(res, vq)
    end

    return res


end