# Zalesak disk test case

using HighOrderReinitialization
using LinearAlgebra
using Plots

grid = Grid((96,96), (0,1))
# sdf for a circle
circle(x) = sqrt((x[1]-0.5)^2 + (x[2]-0.75)^2) - 0.15
# sdf for a single slot in the box
function box(x)
    x = [x[1]-0.5, x[2]-0.65]
    d = abs.(x) - [0.025,0.125]
    return norm([max(d[1],0),max(d[2],0)]) + min(maximum(d),0.0)
end
f(x) = max(circle(x), -box(x)) #circle with a slot
ϕ = set_field(grid,f)
# field_contour(grid,ϕ,levels = [0], clabels=true, show_grid = false,lw = 2)

# velocity of rigid body rotation
u_x(x) = (0.5 - x[2])
u_y(x) = (x[1] - 0.5) 
# convert the velocity to vector field
u_x_v = field_to_vector(set_field(grid,u_x),grid)
u_y_v = field_to_vector(set_field(grid,u_y),grid)

# pre-compute the matices for gradient
Dx, Dy = create_first_derivative_matrix(grid)
ϕ_v = field_to_vector(ϕ, grid)

function evolve(dt, Nt, ϕ_v, u_x_v, u_y_v; reinitialization = true, gap = 2, order = 3, iter = 20, output_figure = true)
    n = 1
    for i = 1:Nt
        # Implicit Euler
        ϕ_v = ( I + dt*(u_x_v.*Dx + u_y_v.*Dy)) \ ϕ_v
        if reinitialization && (i % gap == 0)
            ϕ_v = reinitialization_2D(vector_to_field(ϕ_v, grid), grid, order, iter = iter)
            ϕ_v = field_to_vector(ϕ_v, grid)
        end
        # Output some figures
        if i % 8 == 0 && output_figure
            p = field_contour(grid,vector_to_field(ϕ_v,grid),levels = [0], show_grid=false,clabels=false,lw = 10)
            plot!(p, size=(1000,1000))
            savefig(p, "./save$n")
            n = n + 1
        end
    end
    return vector_to_field(ϕ_v, grid)
end

result = evolve(2*pi/628, 628, ϕ_v, u_x_v, u_y_v, reinitialization = true, gap = 5, order = 3, iter = 300)
p = field_contour(grid,result,levels = [0], clabels=false, show_grid = true,lw = 2)