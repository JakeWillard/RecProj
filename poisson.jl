function stencil1d(m::Int)

    ic = Int(ceil(m/2.0))
    M = zeros((m, m))
    for i=1:m
        for j=1:m
            M[i, j] = (i - ic)^(j-1) / factorial(j-1)
        end
    end
    return inv(M)
end


function periodic_projection(N, m)

    js = [N-m+1:N...]
    append!(js, [1:N...], [1:m...])

    P = sparse([1:N...], [m+1:m+N...], ones(N), N, N+2*m)
    PT = sparse([1:N+2*m...], js, ones(N+2*m), N+2*m, N)

    return P, PT
end


function periodic_derivative_1d(n, m, N)

    Nb = N + 2*m
    D = sparse(zeros(Float64, (Nb, Nb)))
    stencil = stencil1d(m)[n+1,:]
    ic = Int(ceil(m/2.0))

    for i=1:m
        k = i - ic
        vec = stencil[i] * ones(Nb - abs(k))
        D += spdiagm(k => vec)
    end

    P, PT = periodic_projection(N, m)
    return P * D * PT
end


function periodic_derivative(nx, ny, mx, my, Nx, Ny, dx, dy)

    Dx = periodic_derivative_1d(nx, mx, Nx) / dx^nx
    Dy = periodic_derivative_1d(ny, my, Ny) / dy^ny
    return kron(Dy, Dx)
end


function periodic_laplacian(x, y)

    Nx = length(x)
    Ny = length(y)
    dx = x[2] - x[1]
    dy = y[2] - y[1]

    Dxx = periodic_derivative(2, 0, 3, 3, Nx, Ny, dx, dy)
    Dyy = periodic_derivative(0, 2, 3, 3, Nx, Ny, dx, dy)
    return Dxx + Dyy
end


struct PsiSolver

    L :: SparseMatrixCSC{Float64, Int64}
    Dx :: SparseMatrixCSC{Float64, Int64}
    Dy :: SparseMatrixCSC{Float64, Int64}
    Nx :: Int64
    Ny :: Int64

    function PsiSolver(x, y)

        Nx = length(x)
        Ny = length(y)
        dx = x[2] - x[1]
        dy = y[2] - y[1]

        Dx = periodic_derivative(1, 0, 3, 3, Nx, Ny, dx, dy)
        Dy = periodic_derivative(0, 1, 3, 3, Nx, Ny, dx, dy)
        Dxx = periodic_derivative(2, 0, 3, 3, Nx, Ny, dx, dy)
        Dyy = periodic_derivative(0, 2, 3, 3, Nx, Ny, dx, dy)
        L = Dxx + Dyy

        return new(L, Dx, Dy, Nx, Ny)
    end

end


function solve(Bx, By, ps::PsiSolver)

    bx = reshape(Bx, (ps.Nx*ps.Ny, 1))
    by = reshape(By, (ps.Nx*ps.Ny, 1))
    J = ps.Dx*by - ps.Dy*bx

    psi_vec = ps.L \ J
    return reshape(psi_vec, (ps.Nx, ps.Ny))
end



# function psi_periodic(x, y, Bx, By)
#
#     Nx = length(x)
#     Ny = length(y)
#     dx = x[2] - x[1]
#     dy = y[2] - y[1]
#
#     Dx = periodic_derivative(1, 0, 3, 3, Nx, Ny, dx, dy)
#     Dy = periodic_derivative(0, 1, 3, 3, Nx, Ny, dx, dy)
#     Dxx = periodic_derivative(2, 0, 3, 3, Nx, Ny, dx, dy)
#     Dyy = periodic_derivative(0, 2, 3, 3, Nx, Ny, dx, dy)
#     L = lu(Dxx + Dyy)
#
#     bx = reshape(Bx, (Nx*Ny, 1))
#     by = reshape(By, (Nx*Ny, 1))
#     J = Dx*by - Dy*bx
#
#     psi_vec = L \ J
#     psi = reshape(psi_vec, (Nx, Ny))
#
#     return psi
# end
