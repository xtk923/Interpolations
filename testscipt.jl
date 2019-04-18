include("./griddata.jl")


x = -3 .+ 6 .* rand(100);
y = -3 .+ 6 .* rand(100);
v = sin.(x).^4 .* cos.(y);
xq = [-3 .+ rand() .* 6 for i = 1:10];
yq = [-3 .+ rand() .* 6 for i = 1:10];
vq1 = griddata(x, y, v, xq, yq)