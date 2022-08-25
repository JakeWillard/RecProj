using Plots
using HDF5
using Glob
include("./code.jl")


name = "flat"
dat = SimData("./data/$(name)")

anim_psi = animate_flux(dat, frames=10)
anim_E = animate_E(dat, frames=10)
anim_P = animate_P(dat, frames=10)
#anim_rho = animate_rho(dat, frames=10)

gif(anim_psi, "./data/$(name)/$(name)_psi.gif")
gif(anim_E, "./data/$(name)/$(name)_E.gif")
gif(anim_P, "./data/$(name)/$(name)_P.gif")
#gif(anim_rho, "./data/$(name)/$(name)_rho.gif")

