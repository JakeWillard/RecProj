# depends on: Plots, Interpolations

struct FlowField

    f :: Function 
    Nt :: Int64 
    Nl :: Int64

end


function FlowField(Ux, Uy, dx, dy, Nl, Nt; scale=0.01)

    _, _, nt = size(Ux)

    ux = interpolate(Ux, BSpline(Linear()))
    uy = interpolate(Uy, BSpline(Linear()))
    uxe = extrapolate(ux, Flat())
    uye = extrapolate(uy, Flat())
    f(r, t) = begin 
        i = 1 + (r[1] + 0.5)/dx
        j = 1 + (r[2] + 0.5)/dy
        k = 1 + ((t-1)/(Nt-1))*(nt-1)
        scale*Float64[uxe(i,j,k), uye(i,j,k)]
    end

    return FlowField(f, Nt, Nl)
end


function FlowField(dat, Nl, Nt; scale=0.01)

    dx = dat.xvals[2] - dat.xvals[1]
    dy = dat.yvals[2] - dat.yvals[1] 
    return FlowField(dat.Ux, dat.Uy, dx, dy, Nl, Nt, scale=scale)
end



function rk4(x, y, t, flow::FlowField)

    r = Vector([x, y]) 
    k1 = flow.f(r, t)
    k2 = flow.f(r + k1/2, t+0.5)
    k3 = flow.f(r + k2/2, t+0.5)
    k4 = flow.f(r + k3, t+1)
    rnew = r + (k1 + 2*k2 + 2*k3 + k4)/6

    return rnew[1], rnew[2]
end



function streamlines(x0, y0, flow::FlowField)

    Nt = flow.Nl * flow.Ntr

    X = fill(NaN, (flow.Nl, 2, Nt))
    Y = fill(NaN, (flow.Nl, 2, Nt))

    t = 1
    j = 1

    for _=1:flow.Ntr-1
        
        xn, yn = rk4(x0, y0, t, flow)
        X[1,j,t:t+flow.Nl] .= x0 
        X[2,j,t:t+flow.Nl] .= xn
        Y[1,j,t:t+flow.Nl] .= y0 
        Y[2,j,t:t+flow.Nl] .= yn
        t += 1

        for i=3:flow.Nl
            xnew, ynew = rk4(X[i-1,j,t-1], Y[i-1,j,t-1], t-1, flow)
            X[i,j,t:t+flow.Nl] .= xnew
            Y[i,j,t:t+flow.Nl] .= ynew
            t += 1
        end 

        j = (j == 1) ? 2 : 1
    end

    return X, Y
end


function one_streamline(x0, y0, flow::FlowField; gap=0.5)

    X = fill(NaN, (flow.Nl, flow.Nt))
    Y = fill(NaN, (flow.Nl, flow.Nt))
    t = 1

    while t + 2*flow.Nl <= flow.Nt

        xn, yn = rk4(x0, y0, t, flow)
        X[1,t:t+flow.Nl-2] .= x0 
        X[2,t:t+flow.Nl-2] .= xn
        Y[1,t:t+flow.Nl-2] .= y0 
        Y[2,t:t+flow.Nl-2] .= yn
        t += 1

        for i=3:flow.Nl 
            xn, yn = rk4(X[i-1,t-1], Y[i-1,t-1], t-1, flow)
            X[i,t:t+flow.Nl-2] .= xn
            Y[i,t:t+flow.Nl-2] .= yn
            t += 1 
        end 

        t += Int(floor(flow.Nl*gap))
    end

    return X, Y
end


function many_streamlines(Nx, Ny, flow::FlowField; gap=0.5)

    X = fill(NaN, (flow.Nl, Nx*Ny, flow.Nt))
    Y = fill(NaN, (flow.Nl, Nx*Ny, flow.Nt))
    xs = LinRange(-0.5, 0.5, Nx)
    ys = LinRange(-0.5, 0.5, Ny)
    k = 0

    for x in xs 
        for y in ys 
            k += 1 
            Xs, Ys = one_streamline(x, y, flow, gap=gap)
            X[:,k,:] = Xs[:,:]
            Y[:,k,:] = Ys[:,:]
        end 
    end 

    return X, Y 
end

            




function more_streamlines(N, flow::FlowField)

    X = fill(NaN, (flow.Nl, 2*N, flow.Nl*flow.Ntr))
    Y = fill(NaN, (flow.Nl, 2*N, flow.Nl*flow.Ntr))

    for i=1:N
        k = 2*(i-1) + 1
        Xi, Yi = streamlines(flow)
        X[:,k:k+1,:] = Xi[:,:,:]
        Y[:,k:k+1,:] = Yi[:,:,:]
    end

    return X, Y 
end

# TODO: 
# no more randomization: each stream always starts from the same x0 y0. wozniakowski?
# indexing seems a little weird still 
# make the flow field a unit vector? at least fix normalization
#  