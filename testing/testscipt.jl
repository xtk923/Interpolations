include("$(@__DIR__)../src/griddata.jl")
using DelimitedFiles
using Debugger


# SWAP
Mat_X = readdlm("Mat_X.csv")
Mat_Y = readdlm("Mat_Y.csv")
Mat_Z = readdlm("Mat_Z.csv")

idxNan = findall(isnan.(Mat_Z))
idxVal = findall(.!isnan.(Mat_Z))

x = Mat_X[idxVal]
y = Mat_Y[idxVal]
v = Mat_Z[idxVal]
xq = Mat_X[idxNan]
yq = Mat_Y[idxNan]

Mat_Z[idxNan] = griddata(x, y, v, xq, yq)



# Libor
Mat_X = readdlm("Mat_X_Libor.csv")
Mat_Y = readdlm("Mat_Y_Libor.csv")
Mat_Z = readdlm("Mat_Z_Libor.csv");

idxNan = findall(isnan.(Mat_Z));
idxVal = findall(.!isnan.(Mat_Z));

x = Mat_X[idxVal];
y = Mat_Y[idxVal];
v = Mat_Z[idxVal];
xq = Mat_X[idxNan];
yq = Mat_Y[idxNan];

Mat_Z[idxNan] = griddata(x, y, v, xq, yq);