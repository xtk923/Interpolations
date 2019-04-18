include("./griddata.jl")
using .Interpolations_MAT


x = -3 .+ 6 .* rand(1000);
y = -3 .+ 6 .* rand(1000);
v = sin.(x).^4 .* cos.(y);
xq = [-4 .+ rand() .* 8 for i = 1:100];
yq = [-4 .+ rand() .* 8 for i = 1:100];
vq = griddata(x, y, v, xq, yq)