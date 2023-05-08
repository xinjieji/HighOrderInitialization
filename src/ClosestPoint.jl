# This file includes the final step of re-initialization:
# The idea is: for every grid point, first find the closest point in the cloud
# by k-d tree, then use Newton's method to find the exact closest point at the interface

"""
Given a point x, its closest point in the clouds, and the corresponding poly,
using Newton's method to find out the exact closest point at the interface.

- x: the point to be optimized
- x₀: the closest point in the clouds
- poly: the polynomial function of the interface
- r: the radius of the ball
- tol: the tolerance of the Newton's method
- iter: the maximum number of iterations
"""
function newton(x,x₀,poly::Function, r, tol, iter)
    d = length(x)
    λ = (x - x₀)'*ForwardDiff.gradient(poly, x₀)/norm(ForwardDiff.gradient(poly, x₀))
    x_n = x₀
    for k = 1:iter
        g = [x_n - x + λ*ForwardDiff.gradient(poly, x_n); poly(x_n)]
        H = [
            I + λ*ForwardDiff.hessian(poly, x_n)  ForwardDiff.gradient(poly, x_n)
            ForwardDiff.gradient(poly, x_n)'  0
            ]
        if abs(cond(H)) < Inf
            Δ = H\g
            if norm(Δ[1:d]) > 0.5*r
                Δ = 0.5*r*Δ/norm(Δ[1:d])
            end
            x_n1 = x_n - Δ[1:d]
            λ = λ - Δ[d+1]
        else
            # Special situation (H is singular), not tested
            Δ₁ = -(poly(x_n)/norm(ForwardDiff.gradient(poly, x_n))^2)*ForwardDiff.gradient(poly, x_n)
            λ = (x - x_n)'*ForwardDiff.gradient(poly, x_n)/norm(ForwardDiff.gradient(poly, x_n))^2
            Δ₂ = x - x_n - λ*ForwardDiff.gradient(poly, x_n)
            if norm(Δ₂) > 0.1*r
                Δ₂ = 0.1*r*Δ₂/norm(Δ₂)
            end
            x_n1 = x_n + Δ₁ + Δ₂
            print(x₀)
        end
        if norm(x_n1 - x₀) > r
            print("Did not converge with Ball($x₀,$r) for query point $x.")
            return x_n1
        elseif norm(x_n1 - x_n) < tol
            return x_n1
        end
        x_n = x_n1
    end
    return x₀
    # throw("Did not converge in maximum iterations.")
end

"""
Apply the scheme to all the grid points, then generate the new ϕ field.
"""
function closest_point_2D(clouds, Cell_num, C, p, old_ϕ, grid::Grid; iter = 20, r = 2*minimum(grid_spacing(grid)), tol = minimum(grid_spacing(grid))^p)
    kdtree = KDTree(clouds) # build the k-d tree
    ϕ = CircularArray(0.0, grid.dims)
    # loop over every grid point
    for I in eachindex(ϕ)
        x = pos(grid, Tuple(I))
        # find its closest point
        index,dist = knn(kdtree,x,1)
        x₀ = clouds[:,index[1]]
        cell₀ = Int(Cell_num[index[1]])
        # find the corresponding polynomial
        poly(x) = f_inter_2D(x, cell₀, C, p)
        # use Newton's method to find the exact closest point
        x_exact = newton(x,x₀,poly,r,tol,iter)
        # calculate the signed distance 
        ϕ[I] = sign(old_ϕ[I])*norm(x_exact - x)
    end
    return ϕ
end

"""
Apply the scheme to close grid points, then generate the new ϕ field.
"""
function local_closest_point_2D(clouds, cell_index, Cell_num, C, p, old_ϕ, grid::Grid; iter = 20, r = 2*minimum(grid_spacing(grid)), tol = minimum(grid_spacing(grid))^p)
    kdtree = KDTree(clouds) # build the k-d tree
    ϕ = copy(old_ϕ)
    # loop over every point
    for n in 1:length(cell_index)
        for i = -4:5
            for j = -4:5
                I = [cell_index[n][1] + i, cell_index[n][2] + j]
                if I[1] < 1 || I[1] > grid.dims[1] || I[2] < 1 || I[2] > grid.dims[2]
                    continue
                end
                x = pos(grid, Tuple(I))
                # find its closest point
                index,dist = knn(kdtree,x,1)
                x₀ = clouds[:,index[1]]
                cell₀ = Int(Cell_num[index[1]])
                # find the corresponding polynomial
                poly(x) = f_inter_2D(x, cell₀, C, p)
                # use Newton's method to find the exact closest point
                x_exact = newton(x,x₀,poly,r,tol,iter)
                # calculate the signed distance 
                ϕ[I[1],I[2]] = sign(old_ϕ[I[1],I[2]])*norm(x_exact - x)
            end
        end
    end
    return ϕ
end

