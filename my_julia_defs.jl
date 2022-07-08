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

@recipe function f(::Type{Val{:magplot}}, xvals, yvals, Vx, Vy)

    Vmag = sqrt.(Vx.^2 + Vy.^2)
    x := xvals
    y := yvals
    z := Vmag
    seriestype := heatmap
    ()
end












@userplot BPlot
@recipe function f(bp::BPlot)

    dat = bp.args[1]
    t = bp.args[2]

    Bx = dat.Bx[:,:,t]
    By = dat.By[:,:,t]
    B = sqrt.(Bx.^2 + By.^2)
    bx = Bx ./ B
    by = By ./ B
    x = dat.xvals
    y = dat.yvals

    @series begin
        seriestype := heatmap
        x, y, transpose(B)
    end

    Nx = length(x)
    Ny = length(y)
    X = hcat([x for i=1:Ny]...)
    Y = vcat([transpose(y) for i=1:Nx]...)
    xvec = reshape(X, (Nx*Ny,1))
    yvec = reshape(Y, (Nx*Ny, 1))
    bvec = hcat(reshape(bx, (Nx*Ny,1)), reshape(by, (Nx*Ny,1)))

    @series begin
        seriestype := quiver
        xvec, yvec, bvec
    end
end


function plot_B(x, y, Bx, By; scale=0.05, cut=1)

    Nx = length(x)
    Ny = length(y)

    B = sqrt.(Bx.^2 + By.^2)
    bx = Bx ./ B
    by = By ./ B

    h = heatmap(x, y, transpose(B))

    is = 1:cut:Nx
    js = 1:cut:Ny
    Nk = length(is)*length(js)
    println(Nk)
    xs = zeros((2,Nk))
    ys = zeros((2,Nk))

    p = Progress(Nk, 0.5)
    k = 1
    for i=is
        for j=js
            xs[1, k] = x[i]
            xs[2, k] = x[i] + scale*bx[i,j]
            ys[1, k] = y[j]
            ys[2, k] = y[j] + scale*by[i,j]
            k += 1
            next!(p)
        end
    end

    #plot!(h, xs, ys, arrow=true, legend=false)
    return h, xs, js
end



function animate_B(sim::SimData)

    Nt = size(sim.Bx)[3]
    anim = @animate for t=1:Nt
        B = sqrt.(sim.Bx[:,:,t].^2 + sim.By[:,:,t].^2)
        Vx = sim.Ux[:,:,t]
        heatmap(sim.xvals, sim.yvals, transpose(B))#, clims=(-1, 1), aspect_ratio=1)
    end

    return anim
end
