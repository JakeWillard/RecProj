using Plots
using Interpolations
using HDF5
using Glob
using FFTW
include("./my_julia_defs.jl")

dat = SimData("./bigsig")

function plot_stream(x, y, Vx, Vy)

    Nx = length(x)
    Ny = length(y)
    kx = fftfreq(Nx)*Nx
    ky = fftfreq(Ny)*Ny
    kxm = hcat([kx for i=1:Ny]...)
    kym = hcat([ones(Nx)*kyp for kyp in ky]...)
    k2 = kxm .^2 + kym .^2
    k2[1,1] = 1.0

    Vxh = fft(Vx)
    Vyh = fft(Vy)
    Jh = im*(kxm .* Vyh - kym .* Vxh)
    phih = - Jh ./ k2
    phih[1,1] = 0.0

    phi = real.(ifft(phih))
    c = contour(x, y, transpose(phi), color = :black, colorbar=false)
    return c
end
