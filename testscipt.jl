include("./griddata.jl")


x = -3 .+ 6 .* rand(10000);
y = -3 .+ 6 .* rand(10000);
v = sin.(x).^4 .* cos.(y);
xq = [-3 .+ rand() .* 6 for i = 1:100];
yq = [-3 .+ rand() .* 6 for i = 1:100];
vq = griddata(x, y, v, xq, yq)