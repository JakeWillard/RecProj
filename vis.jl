using CairoMakie
using Interpolations
using HDF5
using Glob
include("./my_julia_defs.jl")

dat = SimData("./output")

animate_magnetic_field(dat, "./test.gif", framerate=10)
