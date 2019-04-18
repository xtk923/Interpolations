module Interpolations_MAT
using LinearAlgebra
using VoronoiDelaunay
using DataFrames
using GeometricalPredicates

export griddata
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
    newV, vScale, vShift = normalizedToInterval(v, VoronoiDelaunay.min_coord,
                                                VoronoiDelaunay.max_coord)


    tess = DelaunayTessellation(2*length(x))
    points3D = Point.(newX, newY, newV)
    points = map(p->Point2D(getx(p), gety(p)), points3D)
    push!(tess, points)
    res = []
    df = DataFrame(X = x, Y = y, V = v)
    for i = 1:length(xq)
        normalizedX = xq[i] / xScale - xShift
        normalizedY = yq[i] / yScale - yShift
        p = Point2D(normalizedX, normalizedY)
        theTriangle = locate(tess, p)
        if !isexternal(theTriangle)
            vertices = Array{Float64, 2}(undef, 3,3)
            # index of vertice A in `points`
            idxA = indexin([geta(theTriangle)], points)
            idxB = indexin([getb(theTriangle)], points)
            idxC = indexin([getc(theTriangle)], points)

            vertices[1,:] = [getx(points3D[idxA][1]),
                             gety(points3D[idxA][1]),
                             getz(points3D[idxA][1])]
            vertices[2,:] = [getx(points3D[idxB][1]),
                             gety(points3D[idxB][1]),
                             getz(points3D[idxB][1])]
            vertices[3,:] = [getx(points3D[idxC][1]),
                             gety(points3D[idxC][1]),
                             getz(points3D[idxC][1])]
            vertices[:,1] = reverseNormalizedToInterval(vertices[:,1], xScale,
                                                        xShift)
            vertices[:,2] = reverseNormalizedToInterval(vertices[:,2], yScale,
                                                        yShift)
            vertices[:,3] = reverseNormalizedToInterval(vertices[:,3], vScale,
                                                        vShift)


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
            append!(res, vq)
        else
            append!(res, NaN)
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


"""
        reverse the array of data according to `scale` and `shift`
"""
function reverseNormalizedToInterval(x::Array, scale::Real, shift::Real)

    res = x .+ shift
    res .*= scale

    return res

end

end