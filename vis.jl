using Plots
using Interpolations
using HDF5
using Glob
using FFTW
include("./visualization.jl")

dat = SimData("./output")


function H(x, y, px, py; C1=1.0, C2=1.0)

    phi = C1*sech(sin(y))^2
    chi = C2*sech(sin(y))^2

    return px^2/2 + py^2/(2*(1 + chi)) + phi
end

Ny = 100
Npy = 100
ydist = 1.0
pydist = 10.0

x = 0.0
px = 0.0

y = LinRange(-ydist, ydist, Ny)
py = LinRange(-pydist, pydist, Npy)
ham = hcat([H.(x, y, px, pyp, C1=-10.0, C2=1.0) for pyp in py]...)

contour(y, py, transpose(ham))
