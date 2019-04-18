include("./griddata.jl")
using .Interpolations_MAT


x = -3 .+ 6 .* rand(50000);
y = -3 .+ 6 .* rand(50000);
v = sin.(x).^4 .* cos.(y);
xq = [-3 .+ rand() .* 6 for i = 1:1000];
yq = [-3 .+ rand() .* 6 for i = 1:1000];


vq = griddata(x, y, v, xq, yq)