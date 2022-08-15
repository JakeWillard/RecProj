# depends on: Plots, Interpolations



struct FlowField

    f :: Function 
    Nt :: Int64 
    dt :: Float64

end


function FlowField(Ux, Uy, dx, dy, dt)

    _, _, Nt = size(Ux)
    ux = interpolate(Ux, BSpline(Cubic()))
    uy = interpolate(Uy, BSpline(Cubic()))
    f(r, t) = begin 
        i = 1 + (r[1] + 0.5)/dx
        j = 1 + (r[2] + 0.5)/dy
        k = t/dt 
        Float64[ux(i,j,k), uy(i,j,k)]
    end

    return FlowField(f, Nt, dt)
end



function rk4(x, y, t, flow::FlowField)

    r = Vector([x, y]) 
    k1 = flow.f(r, t)
    k2 = f(r + flow.dt*k1/2, t+flow.dt/2)
    k3 = f(r + flow.dt*k2/2, t+flow.dt/2)
    k4 = f(r + flow.dt*k3, t+flow.dt)
    rnew = r + flow.dt*(k1 + 2*k2 + 2*k3 + k4)/6

    return rnew[1], rnew[2], t + flow.dt
end



function streamline(x0, y0, t0, Nt, Nl, flow::FlowFIeld)

    X = fill(NaN, (Nl, Nt))
    Y = fill(NaN, (Nl, Nt))
    T = fill(NaN, Nt)

    X[1:Nl+1,1] .= x0 
    Y[1:Nl+1,1] .= y0 
    T[1] = t0

    for i=2:Nt 
        xnew, ynew, tnew = rk4(X[i-1,i-1], Y[i-1,i-1], T[i-1], flow)

        X[i:i+Nl,i] .= xnew
        Y[i:i+Nl,i] .= ynew
        T[i] = tnew 
    end

    return X, Y
end 


function animate(Nd, Nt, Nl, flow::FlowField)

    X = fill(NaN, Nl, Nd, flow.Nt)
    Y = fill(NaN, Nl, Nd, flow.Nt)

end






    


    









        

end
    




