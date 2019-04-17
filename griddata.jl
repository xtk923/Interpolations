using LinearAlgebra
using DataFrames
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
        dists = ((xq[i] .- x).^2 + (yq[i] .- y).^2).^0.5
        point = DataFrame(X = xq[i], Y = yq[i], V = [missing], Dist = [0])
        points = DataFrame(X = x, Y = y, V = v, Dist = dists)
        sort!(points, :Dist)
        # find the nearest three by (x,y)

        # if all nearest points are on the left of quest point,
        # drop the left most one.
        vertices = points[1:3, :]
        # the row at which we read the next point from points
        splitRow = 4
        while splitRow <= size(points, 1)
            if !isInTriangle(point, vertices)
                vertices = vertices[1:2, :]
                append!(vertices, points[splitRow, :])
                splitRow += 1
            else
                break
            end
        end
        @debug splitRow
        if splitRow <= size(points, 1) + 1
            vertices = convert(Matrix, vertices[[:X, :Y, :V]])

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
        else
            append!(res, NaN)
        end
    end

    return res


end

"""
        calculate the area of a triangle from its three vertices
        """
function areaTriangle(vertices::DataFrame)
    if size(vertices, 1) != 3
        error("Triangle should have three vertices")
    end

    area = abs(1/2 * (
        vertices.X[1] * (vertices.Y[2] - vertices.Y[3]) +
        vertices.X[2] * (vertices.Y[3] - vertices.Y[1]) +
        vertices.X[3] * (vertices.Y[1] - vertices.Y[2])
    ))
    return area


end


"""
        check if a point is within a triangle
        """
function isInTriangle(point::DataFrame, vertices::DataFrame)
    if size(vertices, 1) != 3
        error("Triangle should have three vertices")
    end

    TA = vcat(vertices[1:2, :], point)
    TB = vcat(vertices[[1,3], :], point)
    TC = vcat(vertices[2:3, :], point)

    AreaA = areaTriangle(TA)
    AreaB = areaTriangle(TB)
    AreaC = areaTriangle(TC)

    AreaOrigin = areaTriangle(vertices)

    return AreaOrigin == AreaA + AreaB + AreaC
end