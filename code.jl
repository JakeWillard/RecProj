#dependencies:
# Plots, HDF5, Glob 


# data structure for simulation output 
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


function flux_function(dat::SimData)

    Mx = zeros(dat.Nx, dat.Ny, dat.Nt)
    My = zeros(dat.Nx, dat.Ny, dat.Nt)

    half_dx = (dat.xvals[2] - dat.xvals[1])/2.0
    half_dy = (dat.yvals[2] - dat.yvals[1])/2.0

    for i=2:dat.Nx
        Mx[i,:,:] = Mx[i-1,:,:] + half_dx*(dat.By[i-1,:,:] + dat.By[i,:,:])
    end

    for j=2:dat.Ny
        My[:,j,:] = My[:,j-1,:] + half_dy*(dat.Bx[:,j-1,:] + dat.Bx[:,j,:])
    end

    return Mx - My 
end


# animate the out of plane vector potential
function animate_flux(dat::SimData; frames=-1, yrange=(-0.5, 0.5))

    psi = flux_function(dat)
    dy = dat.yvals[2] - dat.yvals[1]

    # get index range for range in y
    jmin = 1 + Int(floor((yrange[1] + 0.5)/dy))
    jmax = 1 + Int(floor((yrange[2] + 0.5)/dy))
    yvals = dat.yvals[jmin:jmax]

    # psi_min = minimum(psi[:,jmin:jmax,:])
    # psi_max = maximum(psi[:,jmin:jmax,:])
    # clims = (psi_min, psi_max)

    Nf = (frames == -1) ? dat.Nt : frames

    anim = @animate for t=1:Nf
        contour(dat.xvals, yvals, transpose(psi[:,jmin:jmax,t]), title="psi || timestep: $(t)", color=:black)
    end

    return anim

end


# animating the out of plane E field 
function animate_E(dat::SimData; frames=-1, clims=(-1,1))

    Nf = (frames == -1) ? dat.Nt : frames

    E = dat.Ux .* dat.By - dat.Uy .* dat.Bx
    anim = @animate for t=1:Nf
        heatmap(dat.xvals, dat.yvals, transpose(E[:,:,t]), title="Ez || timestep: $(t)", clims=clims)
    end

    return anim

end


function animate_P(dat::SimData; frames=-1, clims=(0, 1))

    Nf = (frames == -1) ? dat.Nt : frames

    anim = @animate for t=1:Nf
        heatmap(dat.xvals, dat.yvals, transpose(dat.P[:,:,t]), title="P || timestep: $(t)", clims=clims)
    end

    return anim 
end 


function animate_rho(dat::SimData; frames=-1, clims=(0, 1))

    Nf = (frames == -1) ? dat.Nt : frames

    anim = @animate for t=1:Nf
        heatmap(dat.xvals, dat.yvals, transpose(dat.rho[:,:,t]), title="P || timestep: $(t)", clims=clims)
    end

    return anim 
end 
