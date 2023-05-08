# This file includes the functions needed to sample the interpolated interface.
# Strategy: 
# Subdivide each grid cell containing the interface into a 2 × 2 subgrid (2D),
# or a 2 × 2 × 2 subgrid (3D), and in each subgrid, place a point at the center as a seed.
# Then project each seed onto the interface using Newton’s method.

"""
Given a cell, return the 4 seeds in 2D

[seed_1 seed_2 seed_3 seed_4] = seeds_2D(cell_number, cell_index, grid::Grid)
"""
function seeds_2D(cell_num, cell_index, grid::Grid)
    x_ll = pos(grid, cell_index[cell_num]) # lower left corner of the cell
    h = grid_spacing(grid)
    seeds = zeros(2,4)
    seeds[:,1] = x_ll + [h[1]/4, h[2]/4]
    seeds[:,2] = x_ll + [h[1]/4, 3*h[2]/4]
    seeds[:,3] = x_ll + [3*h[1]/4, 3*h[2]/4]
    seeds[:,4] = x_ll + [3*h[1]/4, h[2]/4]
    return seeds
end

"""
Given a cell and its interpolation coefficients, 
use the iterative method to project its seeds to the interface.
- iter: the number of iterations.
- tol: ||x_{n+1} - x_n|| < tol * (h_min/2), 0.01 is ok by default.
"""
function iterate_proj_2D(cell_num, cell_index, C, p, grid::Grid; iter = 10, tol = 1e-2)
    seeds = seeds_2D(cell_num, cell_index, grid)
    f(x) = f_inter_2D(x, cell_num, C, p)
    d,ns = size(seeds)
    keep = []
    for n in 1:ns
        seeds_0 = copy(seeds[:,n])
        for i in 1:iter
            tmp = ForwardDiff.gradient(f, seeds[:,n]) * f(seeds[:,n]) / norm(ForwardDiff.gradient(f, seeds[:,n]))^2
            if norm(tmp) < tol*minimum(grid_spacing(grid)/2)
                break
            end
            seeds[:,n] -= tmp
        end
        # Drop the seed if it is far away from the original position
        if norm(seeds[:,n] - seeds_0) < minimum(grid_spacing(grid)/2)
            push!(keep, n)
        end
    end
    return seeds[:, keep]
end

"""
Sample all the points given the cell_index, interpolation coefficients, and the grid.

clouds = sample_2D(cell_index, C, p, grid::Grid; iter = 10, tol = 1e-2)
"""
function sample_2D(cell_index, C, p, grid::Grid; iter = 10, tol = 1e-2)
    Seeds = []
    Cell_num = []
    for i in 1:length(cell_index)
        seeds = iterate_proj_2D(i, cell_index, C, p, grid, iter = iter, tol = tol)
        cell_num = zeros(1,size(seeds,2)) .+ i
        push!(Seeds, seeds)
        push!(Cell_num, cell_num)
    end
    clouds = x_wrapper(Seeds)
    Cell_num = x_wrapper(Cell_num)
    return clouds, Cell_num
end

