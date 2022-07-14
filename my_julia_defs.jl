struct SimData

    Bx :: Array{Float32, 3}
    By :: Array{Float32, 3}
    rho :: Array{Float32, 3}
    P :: Array{Float32, 3}
    Ux :: Array{Float32, 3}
    Uy :: Array{Float32, 3}

    xvals :: Vector{Float32}
    yvals :: Vector{Float32}

    Nx :: Int32
    Ny :: Int32
    Nt :: Int32

    function SimData(datafolder)

        files = glob("$(datafolder)/*.athdf")

        fid = h5open(files[1], "r")
        xvals = Vector(fid["x1v"][:,1])
        yvals = Vector(fid["x2v"][:,1])
        Nx = length(xvals)
        Ny = length(yvals)
        Nt = length(files)
        close(fid)

        Bx = zeros(Float32, (Nx, Ny, Nt))
        By = zeros(Float32, (Nx, Ny, Nt))
        rho = zeros(Float32, (Nx, Ny, Nt))
        P = zeros(Float32, (Nx, Ny, Nt))
        Ux = zeros(Float32, (Nx, Ny, Nt))
        Uy = zeros(Float32, (Nx, Ny, Nt))

        for t=1:Nt
            fid = h5open(files[t], "r")
            Bx[:,:,t] = fid["B"][:,:,1,1,1]
            By[:,:,t] = fid["B"][:,:,1,1,2]
            rho[:,:,t] = fid["prim"][:,:,1,1,1]
            P[:,:,t] = fid["prim"][:,:,1,1,2]
            Ux[:,:,t] = fid["prim"][:,:,1,1,3]
            Uy[:,:,t]=  fid["prim"][:,:,1,1,4]
            close(fid)
        end

        new(Bx, By, rho, P, Ux, Uy, xvals, yvals, Nx, Ny, Nt)
    end

end

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
