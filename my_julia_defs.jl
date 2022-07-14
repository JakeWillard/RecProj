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
