using Plots
using Interpolations
using HDF5
using Glob
using FFTW
using PlutoUI
include("./visualization.jl")

dat = SimData("./output")
