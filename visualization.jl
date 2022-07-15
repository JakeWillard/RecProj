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


function stream_function(Vx::Matrix{Float32}, Vy::Matrix{Float32})

    Nx, Ny = size(Vx)
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

    return real.(ifft(phih))
end


function stream_function(Vx::Array{Float32, 3}, Vy::Array{Float32, 3})

    Nx, Ny, Nt = size(Vx)
    kx = fftfreq(Nx)*Nx
    ky = fftfreq(Ny)*Ny
    kxm = hcat([kx for i=1:Ny]...)
    kym = hcat([ones(Nx)*kyp for kyp in ky]...)
    k2 = kxm .^2 + kym .^2
    k2[1,1] = 1.0

    phi = zeros((Nx, Ny, Nt))
    phih = zeros(ComplexF32, (Nx, Ny))

    # compute initial phi
    FFT = plan_fft(Vx[:,:,1])
    Jh = im*(kxm .* (FFT * Vy[:,:,1]) - kym .* (FFT * Vx[:,:,1]))
    phih[:,:] = - Jh ./ k2
    phih[1,1] = 0.0
    IFFT = plan_ifft(phih)
    phi[:,:,1] = real.(IFFT * phih)

    for t=2:Nt

        Jh = im*(kxm .* (FFT * Vy[:,:,t]) - kym .* (FFT * Vx[:,:,t]))
        phih[:,:] = -Jh ./ k2
        phih[1,1] = 0.0
        phi[:,:,t] = real.(IFFT * phih)

    end

    return phi
end
