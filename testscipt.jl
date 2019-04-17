using LinearAlgebra
using VoronoiDelaunay


x = -3 .+ 6 .* rand(50);
y = -3 .+ 6 .* rand(50);
v = sin.(x).^4 .* cos.(y);
xq = [-3 .+ rand() .* 6 for i = 1:10]
yq = [-3 .+ rand() .* 6 for i = 1:10]

max = max_coord
min = min_coord

newX, xScale, xShift = normalizedToInterval(x, min, max)
newY, yScale, yShift = normalizedToInterval(y, min, max)

backX = reverseNormalizedToInterval(newX, xScale, xShift)
backX == x

tess = DelaunayTessellation(50)

points = Point.(x, y)

push!(tess, points)

triangles = []
for dTriangle in tess
    push!(triangles, dTriangle)
end

function normalizedToInterval(x::Array, min::Real, max::Real)
    width = max - min
    dataMax = maximum(x)
    dataMin = minimum(x)
    scale  = (dataMax- dataMin) / width
    res = x ./ scale
    shift = minimum(x) - min
    res .-= shift
    return x, scale, shift
end


function reverseNormalizedToInterval(x::Array, scale::Real, shift::Real)

    res = x .+ shift
    res .*= scale

    return x

end
