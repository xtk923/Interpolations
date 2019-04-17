module Interpolations_MAT
using LinearAlgebra
using DataFrames
using VoronoiDelaunay
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

    tess = DelaunayTessellation(2*length(x))
    points = Point.(newX, newY)
    push!(tess, points)
    # tess is a Delaunay Trianglation
    triangles = []
    for t in tess
        push!(triangles, Primitive(geta(t), getb(t), getc(t)))
    end

    res = []
    for i = 1:length(xq)
        normalizedX = xq[i] / xScale - xShift
        normalizedY = yq[i] / yScale - yShift
        p = Point2D(normalizedX, normalizedY)
        theTriangle = nothing
        for t in triangles
            if isInTriangle(p, t)
                theTriangle = t
                break
            end
        end

        if theTriangle != nothing
            skipping = false
            vertices = Array{Float64, 2}(undef, 3,3)
            vertices[1,1] = getx(geta(theTriangle))
            vertices[2,1] = getx(getb(theTriangle))
            vertices[3,1] = getx(getc(theTriangle))
            vertices[1,2] = gety(geta(theTriangle))
            vertices[2,2] = gety(getb(theTriangle))
            vertices[3,2] = gety(getc(theTriangle))
            vertices[:,1] = reverseNormalizedToInterval(vertices[:,1], xScale, xShift)
            vertices[:,2] = reverseNormalizedToInterval(vertices[:,2], yScale, yShift)

            for r = 1:3
                idxAccX = findall(x .≈ vertices[r, 1])
                idxAccY = findall(y .≈ vertices[r, 2]) 
                idxForZ = intersect(idxAccX, idxAccY)
                if length(idxForZ) == 1
                    vertices[r, 3] = v[idxAccX[1]]
                else
                    skipping = true
                end
            end

            if !skipping
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
        else
            append!(res, NaN)
        end
    end

    return res

end

"""
        calculate the area of a triangle from its three vertices
"""
function areaTriangle(t::GeometricalPredicates.UnOrientedTriangle{Point2D})
    area = abs(1/2 * (
        getx(geta(t)) * (gety(getb(t)) - gety(getc(t))) +
        getx(getb(t)) * (gety(getc(t)) - gety(geta(t))) +
        getx(getc(t)) * (gety(geta(t)) - gety(getb(t)))
    ))
    return area
end


"""
        check if a point is within a triangle
    """
function isInTriangle(p::Point2D,
                      t::GeometricalPredicates.UnOrientedTriangle{Point2D})
    TA = Primitive(p, getb(t), getc(t))

    TB = Primitive(p, geta(t), getc(t))

    TC = Primitive(p, geta(t), getb(t))

    AreaA = areaTriangle(TA)
    AreaB = areaTriangle(TB)
    AreaC = areaTriangle(TC)

    AreaOrigin = areaTriangle(t)

    return AreaOrigin ≈ AreaA + AreaB + AreaC
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