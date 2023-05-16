# Test the static re-initialization method for an ellipse interface
using HighOrderReinitialization
using LinearAlgebra

f(x) = (1-exp(-(x[1]-0.3)^2 - (x[2]-0.3)^2))*(sqrt(4*x[1]^2+9*x[2]^2)-1)

function exact_ellipse(x, p)
    x = abs.(x)
    if x[1]>x[2]
        x = reverse(x)
        p = reverse(p)
    end
    l = p[2]^2 - p[1]^2
    m = p[1]*x[1]/l
    m2 = m^2
    n = p[2]*x[2]/l
    n2 = n^2
    c = (m2 + n2 - 1)/3
    c3 = c^3
    q = c3 + m2*n2*2
    d = c3 + m2*n2
    g = m + m*n2
    if d<0
        h = acos(q/c3)/3
        s = cos(h)
        t = sin(h)*sqrt(3)
        rx = sqrt(-c*(s + t + 2) + m2)
        ry = sqrt(-c*(s - t + 2) + m2)
        co = (ry + sign(l)*rx + abs(g)/(rx*ry) - m)/2
    else
        h = 2*m*n*sqrt(d)
        s = sign(q + h)*abs(q + h)^(1/3)
        u = sign(q - h)*abs(q - h)^(1/3)
        rx = -s - u - c*4 + 2*m2
        ry = (s - u)*sqrt(3)
        rm = sqrt(rx^2 + ry^2)
        co = (ry/sqrt(rm-rx) + 2*g/rm-m)/2
    end
    r = p.*[co, sqrt(1-co^2)]
    if x == [0,0]
        return -minimum(p)
    else
        return norm(r-x)*sign(x[2]-r[2])
    end
end

function error(Nx::Vector, p; iter = 20, tol = 1e-13)
    err_2 = Float64[]
    err_inf = Float64[]
    for i in 1:length(Nx)
        grid = Grid((Nx[i],Nx[i]), (-3/4,3/4))
        ϕ = set_field(grid,f)
        ϕ_new = reinitialization_2D(ϕ, grid, p, iter = iter, tol=tol)
        ϕ_e = set_field(grid,(x)->exact_ellipse(x,[1/2,1/3]))
        push!(err_2, norm(ϕ_new - ϕ_e)/Nx[i])
        push!(err_inf, maximum(abs, ϕ_new - ϕ_e))
    end
    return [Nx err_2], [Nx err_inf]
end

function order(err)
    N = err[:,1]
    err_2 = err[:,2]
    order_2 = log.(err_2[1:end-1]./err_2[2:end])./log.(N[2:end]./N[1:end-1])
    return order_2
end

Nx = [48, 96, 192, 384]
err_2, err_inf = error(Nx, 4, iter=20, tol =1e-13)
order_2, order_inf = order(err_2), order(err_inf)