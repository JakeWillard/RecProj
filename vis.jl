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

    return (1 - C1*phi)*sqrt(1 + px^2 + (1 + chi)*py^2)
end

Ny = 100
Npy = 100
ydist = 1.0
pydist = 1.0

x = 0.0
px = 0.0

y = LinRange(-ydist, ydist, Ny)
py = LinRange(-pydist, pydist, Npy)
ham = hcat([H.(x, y, px, pyp, C1=-0.9, C2=10.0) for pyp in py]...)

contour(y, py, transpose(ham))
