using CairoMakie
using Interpolations
using HDF5
using Glob
include("./my_julia_defs.jl")

dat = SimData("./output")

animate_magnetic_field(dat, "./test.gif")

fig, ax, hm = heatmap(dat.xvals, dat.yvals, dat.rho[:,:,10])
heatmap!(ax, dat.xvals, dat.yvals, dat.rho[:,:,20])

hm = heatmap(dat.xvals, dat.yvals, dat.rho[:,:,40])
