using CairoMakie
using Interpolations
using HDF5
using Glob
include("./my_julia_defs.jl")

dat = SimData("./output")
