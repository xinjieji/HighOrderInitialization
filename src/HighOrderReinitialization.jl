module HighOrderReinitialization

using Meshing
using GeometryBasics
using ImplicitPlots
using Plots
using CircularArrays
using LinearAlgebra
using ForwardDiff
using NearestNeighbors
using OffsetArrays
using SparseArrays

include("Grids.jl")
export Grid, grid_spacing, pos, set_field, field_to_vector, find_index, vector_to_field

include("Plots.jl")
export implicit_plot_3D, level_set_contour, field_contour, implicit_plot_2D

include("PiecewiseApprox.jl")
export x_wrapper, cal_coeff_2D, f_inter_2D

include("InterfaceSample.jl")
export sample_2D

include("ClosestPoint.jl")
export closest_point_2D, local_closest_point_2D

include("Evlove.jl")
export create_first_derivative_matrix

"""
General function for the reinitialization
- ϕ: the field to be reinitialized
- grid: the grid
- p: the order of the polynomial
"""
function reinitialization_2D(ϕ, grid::Grid, p; iter = 20, r = 1.5*minimum(grid_spacing(grid)), tol = minimum(grid_spacing(grid))^p)
    # Step 1: interpolate the cells that cover the interface
    cell_index, X, C = cal_coeff_2D(grid,ϕ,p)
    # Step 2: sample the interface
    clouds, Cell_num = sample_2D(cell_index, C, p, grid)
    # Step 3: kd-tree to find the nearest points
    ϕ_new = closest_point_2D(clouds, Cell_num, C, p, ϕ, grid, iter = iter, r = r, tol = tol)
    # ϕ_new = local_closest_point_2D(clouds, cell_index, Cell_num, C, p, ϕ, grid, iter = iter, r = r, tol = tol)
    return ϕ_new
end
export reinitialization_2D

end # module HighOrderInitialization
