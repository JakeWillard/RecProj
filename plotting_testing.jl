using CairoMakie
using Interpolations
using HDF5
using Glob
include("./my_julia_defs.jl")

function plot_magnetic_field(x, y, Bx, By)

    Bx_interp = interpolate(Bx, BSpline(Linear()))
    By_interp = interpolate(By, BSpline(Linear()))
    b(p::Point2) = begin
        i = 1 + (p[1] - x[1]) / (x[2] - x[1])
        j = 1 + (p[2] - y[1]) / (y[2] - y[1])
        bx = Bx_interp(i,j)
        by = By_interp(i,j)
        Point2(bx, by)
    end

    streamplot(b, x[1]..x[end], y[1]..y[end])
end

plot_magnetic_field(x, y, Bx, By)




# function magplot(x, y, Bx, By, ps::PsiSolver)
#
#     B = sqrt.(Bx.^2 + By.^2)
#     psi = solve(Bx, By, ps)
#
#     h = heatmap(x, y, transpose(B))
#     contour!(h, x, y, transpose(psi), color=:black)
# end


# @userplot MagPlot
# @recipe function f(mp::MagPlot)
#
#     x, y, Bx, By = mp.args
#     B = sqrt.(Bx.^2 + By.^2)
#
#     psi = psi_periodic(x, y, Bx, By)
#
#     # @series begin
#     #     seriestype := heatmap
#     #     x, y, transpose(B)
#     # end
#
#     @series begin
#         seriestype = :contour
#         #color := :black
#         x, y, transpose(psi)
#     end
#
# end




#
#
# @userplot VectorPlot
# @recipe function f(vp::VectorPlot, Nax=20, Nay=20, alen=0.05)
#
#     x, y, Vx, Vy = vp.args
#     Nx = length(x)
#     Ny = length(y)
#
#     V = sqrt.(Vx.^2 + Vy.^2)
#     vx = Vx ./ V
#     vy = Vy ./ V
#
#     cutx = Int64(floor(Nx / Nax))
#     cuty = Int64(floor(Ny / Nay))
#     is = Int64[1:cutx:Nx...]
#     js = Int64[1:cuty:Ny...]
#     Nk = length(is)*length(js)
#
#     X = zeros((2, Nk))
#     Y = zeros((2, Nk))
#
#     k = 1
#     for i in is
#         for j in js
#             X[1,k] = x[i]
#             X[2,k] = x[i] + alen*vx[i,j]
#             Y[1,k] = y[j]
#             Y[2,k] = y[j] + alen*vy[i,j]
#             k += 1
#         end
#     end
#
#     @series begin
#         seriestype := heatmap
#         aspect_ratio --> 1
#         x, y, transpose(V)
#     end


    # @series begin
    #     seriestype := :line
    #     legend := false
    #     color := :black
    #     arrow := true
    #     aspect_ratio --> 1
    #     X, Y
    # end

# end





#
# @recipe function f(::Type{Val{:magplot}}, xvals, yvals, Vx, Vy)
#
#     Vmag = sqrt.(Vx.^2 + Vy.^2)
#     x := xvals
#     y := yvals
#     z := Vmag
#     seriestype := heatmap
#     ()
# end
#
# @userplot BPlot
# @recipe function f(bp::BPlot)
#
#     dat = bp.args[1]
#     t = bp.args[2]
#
#     Bx = dat.Bx[:,:,t]
#     By = dat.By[:,:,t]
#     B = sqrt.(Bx.^2 + By.^2)
#     bx = Bx ./ B
#     by = By ./ B
#     x = dat.xvals
#     y = dat.yvals
#
#     @series begin
#         seriestype := heatmap
#         x, y, transpose(B)
#     end
#
#     Nx = length(x)
#     Ny = length(y)
#     X = hcat([x for i=1:Ny]...)
#     Y = vcat([transpose(y) for i=1:Nx]...)
#     xvec = reshape(X, (Nx*Ny,1))
#     yvec = reshape(Y, (Nx*Ny, 1))
#     bvec = hcat(reshape(bx, (Nx*Ny,1)), reshape(by, (Nx*Ny,1)))
#
#     @series begin
#         seriestype := quiver
#         xvec, yvec, bvec
#     end
# end
#
#
# function plot_B(x, y, Bx, By; scale=0.05, cut=1)
#
#     Nx = length(x)
#     Ny = length(y)
#
#     B = sqrt.(Bx.^2 + By.^2)
#     bx = Bx ./ B
#     by = By ./ B
#
#     h = heatmap(x, y, transpose(B))
#
#     is = 1:cut:Nx
#     js = 1:cut:Ny
#     Nk = length(is)*length(js)
#     println(Nk)
#     xs = zeros((2,Nk))
#     ys = zeros((2,Nk))
#
#     p = Progress(Nk, 0.5)
#     k = 1
#     for i=is
#         for j=js
#             xs[1, k] = x[i]
#             xs[2, k] = x[i] + scale*bx[i,j]
#             ys[1, k] = y[j]
#             ys[2, k] = y[j] + scale*by[i,j]
#             k += 1
#             next!(p)
#         end
#     end
#
#     #plot!(h, xs, ys, arrow=true, legend=false)
#     return h, xs, js
# end
#
#
#
# function animate_B(sim::SimData)
#
#     Nt = size(sim.Bx)[3]
#     anim = @animate for t=1:Nt
#         B = sqrt.(sim.Bx[:,:,t].^2 + sim.By[:,:,t].^2)
#         Vx = sim.Ux[:,:,t]
#         heatmap(sim.xvals, sim.yvals, transpose(B))#, clims=(-1, 1), aspect_ratio=1)
#     end
#
#     return anim
# end



function animate_magnetic_field(dat::SimData, path; framerate=30)

    x = dat.xvals
    y = dat.yvals
    dx = x[2] - x[1]
    dy = y[2] - y[1]

    fig = Figure()
    ax1 = Axis(fig[1, 1:2], title="Magnetic Field", xlabel="X", ylabel="Y")
    ax2 = Axis(fig[2, 1:2], title="Velocity", xlabel="X", ylabel="Y")
    ax3 = Axis(fig[3, 1], title="Density", xlabel="X", ylabel="Y")
    ax4 = Axis(fig[3, 2], title="Pressure", xlabel="X", ylabel="Y")

    record(fig, path, [1:dat.Nt...], framerate=framerate) do t

        Bx = interpolate(dat.Bx[:,:,t], BSpline(Linear()))
        By = interpolate(dat.By[:,:,t], BSpline(Linear()))
        Ux = interpolate(dat.Ux[:,:,t], BSpline(Linear()))
        Uy = interpolate(dat.Uy[:,:,t], BSpline(Linear()))
        rho = dat.rho[:,:,t]
        P = dat.P[:,:,t]

        b(p::Point2) = begin
            i = 1 + (p[1] - x[1]) / dx
            j = 1 + (p[2] - y[1]) / dy
            bx = Bx(i,j)
            by = By(i,j)
            Point2(bx, by)
        end

        vel(p::Point2) = begin
            i = 1 + (p[1] - x[1]) / dx
            j = 1 + (p[2] - y[1]) / dy
            velx = Ux(i,j)
            vely = Uy(i,j)
            Point2(1.0, 1.0)
        end


        empty!(ax1)
        empty!(ax2)
        empty!(ax3)
        empty!(ax4)

        streamplot!(ax1, b, x[1]..x[end], y[1]..y[end])
        streamplot!(ax2, vel, x[1]..x[end], y[1]..y[end])
        heatmap!(ax3, x, y, rho)
        heatmap!(ax4, x, y, P)


    end
end
