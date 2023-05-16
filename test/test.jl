using HighOrderInitialization

grid = Grid((24,24), (-3/4,3/4))
#f(x) = (x[1]^4 + x[2]^4 - 1) * (x[1]^2 + x[2]^2 - 2) + x[1]^5 * x[2]
#f(x) = 4*x[1]^2 + 9*x[2]^2 - 1
f(x) = (1-exp(-(x[1]-0.3)^2 - (x[2]-0.3)^2))*(sqrt(4*x[1]^2+9*x[2]^2)-1)
ϕ = set_field(grid,f)

# Could plot the figure to check the level set function:
 levels = [-0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3]
# field_contour(grid,ϕ,levels=levels, clabels=false)
level_set_contour(grid,f,levels=levels, clabels=false)




"""
Main steps:
"""

# Step 1: interpolate the cells that cover the interface
cell_index, X, C = cal_coeff_2D(grid,ϕ,3)

# Could plot the figure to check the interpolation points:
# data = x_wrapper(X)
# implicit_plot_2D(grid,f, dataset=data, markersize = 2)
# level_set_contour(grid,f)

# Step 2: sample the interface
clouds, Cell_num = sample_2D(cell_index, C, 3, grid)

# Could plot the figure to check the sampled points:
# implicit_plot_2D(grid,f, dataset=clouds, markersize = 2)

# Step 3: kd-tree and Newton's method to find the nearest points
ϕ_new = closest_point_2D(clouds, Cell_num, C, 3, ϕ, grid)

# Could plot the signed distance function:
# field_contour(grid,ϕ_new,levels=levels, clabels=true)