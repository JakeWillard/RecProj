using Plots
using RecipesBase
using HDF5
using Glob
using ProgressMeter
include("./my_julia_defs.jl")

dat = SimData("./bigsig")
gif(animate_B(dat), "./test.gif", fps=10)
plot(dat.xvals, dat.yvals, dat.Bx[:,:,1], dat.By[:,:,1], t=:magplot)
