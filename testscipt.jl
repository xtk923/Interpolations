include("./griddata.jl")


x = -3 .+ 6 .* rand(100);
y = -3 .+ 6 .* rand(100);
v = sin.(x).^4 .* cos.(y);
xq = [-3 .+ rand() .* 6 for i = 1:10];
yq = [-3 .+ rand() .* 6 for i = 1:10];
vq1 = griddata(x, y, v, xq, yq)


Mat_X = repeat(collect(1:7)',  6)
Mat_Y = repeat(collect(1:6), 1, 7)
Mat_Z = Mat_X .+ Mat_Y

idxNan = collect(1:5:42)
idxVal = filter(x->!(x in idxNan), 1:42)
x = Mat_X[idxVal]
y = Mat_Y[idxVal]
v = Mat_Z[idxVal]
xq = [1,3,4,5]
yq = [5,2.5, 3.5, 1]

vq = griddata(x, y, v, xq, yq)