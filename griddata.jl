using LinearAlgebra
using VoronoiDelaunay
using GeometricalPredicates

"""
        griddata with linear interpolation which works similar to its matlab counterpart

        `griddata(x::Array{Real, 1}, y::Array{Real, 1}, v::Array{Real, 1}, xq::Array{Real, 1}, yq::Array{Real, 1})`
"""
function griddata(x::Array, y::Array, v::Array,
                  xq, yq)
    if !(length(x) == length(y) && length(y) == length(v))
        error("x, y and v must be in same dimension.")
    end
    if !(length(xq) == length(yq))
        error("xq and yq must be in same dimension.")
    end

    newX, xScale, xShift = normalizedToInterval(x, VoronoiDelaunay.min_coord,
                                                VoronoiDelaunay.max_coord)
    newY, yScale, yShift = normalizedToInterval(y, VoronoiDelaunay.min_coord,
                                                VoronoiDelaunay.max_coord)

    tess = DelaunayTessellation(2*length(x))
    points = Point.(newX, newY)
    push!(tess, points)
    dict = Dict(zip(points, collect(1:length(x))))
    res = Array{Float64, 1}(undef, length(xq))
    for i = 1:length(xq)
        normalizedX = xq[i] / xScale - xShift
        normalizedY = yq[i] / yScale - yShift
        p = Point2D(normalizedX, normalizedY)
        theTriangle = locate(tess, p)
        if !isexternal(theTriangle)
            vertices = Array{Float64, 2}(undef, 3,3)
            # index of vertice A in `points`
            idxA = dict[geta(theTriangle)]
            idxB = dict[getb(theTriangle)]
            idxC = dict[getc(theTriangle)]
            vertices[1,:] = [x[idxA], y[idxA], v[idxA]]
            vertices[2,:] = [x[idxB], y[idxB], v[idxB]]
            vertices[3,:] = [x[idxC], y[idxC], v[idxC]]

            # with three nearest points, get the expression of the plane,
            # then find the z value of f(xq[i], yq[i])
            D = det(vertices)
            D = D == 0 ? eps() : D
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
            res[i] = vq
        else
            res[i] = NaN
        end
    end

    return res

end

"""
        Scale and shift an array of data to fit within
        the interval defined by `min` and `max`
"""
function normalizedToInterval(x::Array, min::Real, max::Real)
    width = max - min
    dataMax = maximum(x)
    dataMin = minimum(x)
    scale  = (dataMax- dataMin) / width
    res = x ./ scale
    shift = minimum(res) - min
    res .-= shift
    return res, scale, shift
end
