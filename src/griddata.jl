using LinearAlgebra
using DataFrames
using Triangle
using ProgressMeter

"""
        griddata with linear interpolation which works similar to its matlab counterpart

        `griddata(x::Array{Real, 1}, y::Array{Real, 1}, v::Array{Real, 1}, xq::Array{Real, 1}, yq::Array{Real, 1})`
"""
function griddata(x::Array, y::Array, v::Array, xq, yq)

    if !(length(x) == length(y) && length(y) == length(v))
        error("x, y and v must be in same dimension.")
    end
    if !(length(xq) == length(yq))
        error("xq and yq must be in same dimension.")
    end

    points = [Float64.(x) Float64.(y) Float64.(v)]
    points_map = collect(1:length(x))

    triangles = Triangle.basic_triangulation(points[:, 1:2], points_map)
    res = Array{Float64, 1}(undef, length(xq))
    dists = Array{Float64, 1}(undef, length(x))
    df = DataFrame(idx = points_map, X = x, Y = y, V = v, Dist = dists)
    @showprogress 1 for i = 1:length(xq)
        p = [xq[i]  yq[i]]
        # get the distances
        for i = 1:length(dists)
            dists[i] = ((p[1] - points[i, 1])^2 + (p[2] - points[i, 2])^2)^0.5
        end
        df.Dist = dists
        # @time sort!(df, :Dist);
        filter(row->row[:Dist] < 2, df)
        neighbours = filter(row->row[:Dist] < 2, df).idx
        candidates = []
        for i = 1:length(neighbours)
            candidates = vcat(candidates, filter(t-> neighbours[i] in t, triangles))
        end
        # @time df = sort(df, :idx);
        theTriangle = locate(candidates, p, df)
        if theTriangle != nothing
            vertices = Array{Float64, 2}(undef, 3,3)
            # index of vertice A in `points`
            idxA = theTriangle[1]
            idxB = theTriangle[2]
            idxC = theTriangle[3]
            vertices[1,:] = [points[idxA, 1], points[idxA, 2], points[idxA, 3]]
            vertices[2,:] = [points[idxB, 1], points[idxB, 2], points[idxB, 3]]
            vertices[3,:] = [points[idxC, 1], points[idxC, 2], points[idxC, 3]]
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


function areaTriangle(p1::Array, p2::Array, p3::Array)
    area = abs(1/2 * (
        p1[1] * (p2[2] - p3[2]) +
        p2[1] * (p3[2] - p1[2]) +
        p3[1] * (p1[2] - p2[2])
    ))

    return area
end


function inTriangle(p::Array, t::Array{Int64, 1}, df::DataFrame)
    pA = convert(Matrix, filter(row->row[:idx] == t[1], df)[[:X, :Y]])
    pB = convert(Matrix, filter(row->row[:idx] == t[2], df)[[:X, :Y]])
    pC = convert(Matrix, filter(row->row[:idx] == t[3], df)[[:X, :Y]])

    areaA = areaTriangle(p, pB, pC)
    areaB = areaTriangle(p, pA, pC)
    areaC = areaTriangle(p, pA, pB)
    areaOriginal = areaTriangle(pA, pB, pC)
    return areaOriginal ≈ areaA + areaB + areaC
end

function locate(candidates::Array, p::Array, df::DataFrame)

    for t in candidates
        if inTriangle(p, t, df)
            return t
        end
    end

    return nothing
end


