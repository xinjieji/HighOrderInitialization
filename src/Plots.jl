import GLMakie: Axis3, Figure, mesh!, activate!, wireframe!, scale!, xlims!, ylims!, zlims!, scatter!, scatter

"""
Plot implicitly defined 3D surface, scatter points can be added.

Example of usage:
> f(x) = x[1]^2+x[2]^2-x[3]^2-1
> fig = implicit_plot_3D(f, xlims=(-3,3), ylims=(-3,3), zlims=(-3,3))

It can simoultaneously plot a dataset of points [x y z]:
> fig = implicit_plot_3D(f, dataset = [x y z], markersize = 5, xlims=(-3,3), ylims=(-3,3), zlims=(-3,3))
"""
function implicit_plot_3D(
    f;
    dataset = [0.0 0.0 0.0],
    markersize = 0,
    pointColor = :steelblue,
    x_min = -3.0,
    xmin = x_min,
    x_max = 3.0,
    xmax = x_max,
    y_min = xmin,
    ymin = y_min,
    y_max = xmax,
    ymax = y_max,
    z_min = xmin,
    zmin = z_min,
    z_max = xmax,
    zmax = z_max,
    xlims = (xmin, xmax),
    ylims = (ymin, ymax),
    zlims = (zmin, zmax),
    surfaceColor = :grey,
    transparency = true,
    samples=(35,35,35),
    shading = true,
    wireframe=false,
    MarchingModeIsCubes=true,
    zcolormap=nothing,
    kwargs...
)
    activate!()

    try
        f([1,1,1])
    catch
        throw("The specified function does not take 3 arguments.")
    end

    fig = Figure()
    Axis3(fig[1,1])
    implicit_mesh = Mesh(f,Rect(Vec(xlims[1], ylims[1],zlims[1]),
                            Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingModeIsCubes ? MarchingCubes() : MarchingTetrahedra())
    if wireframe
        wireframe!(implicit_mesh, shading=shading, color=surfaceColor, transparency=transparency)
    else
        vertices = decompose(Point{3, Float64}, implicit_mesh)
        triangles = decompose(TriangleFace{Int}, implicit_mesh)
        if isnothing(zcolormap)
            mesh!(vertices, triangles, shading=shading, color=surfaceColor, transparency=transparency, kwargs...)
        else
            colors = [v[3] for v in vertices]
            mesh!(vertices, triangles, shading=shading, color=surfaceColor, transparency=transparency, colormap=zcolormap, kwargs...)
        end        
    end
    scatter!(dataset, color=pointColor, markersize=markersize)

    return(fig)
end

"""
Plot a 2D level-set countour
level_set_contour(grid, func)
"""
function level_set_contour(grid::Grid, func; show_grid=true, 
    levels=6,
    color=:turbo,
    clabels=false,
    cbar=false,
    lw=2,
    legend=false,
    aspect_ratio=:equal,
    axis = false,
    plot_grid = false,
    kwargs...
    )
    if length(grid.dims) != 2
        throw("The grid must be 2D.")
    end

    newgrid = Grid((200,200), grid.size)
    x = grid.size[1]:grid_spacing(grid)[1]:grid.size[2]-grid_spacing(grid)[1]
    y = grid.size[1]:grid_spacing(grid)[2]:grid.size[2]-grid_spacing(grid)[2]
    newx = newgrid.size[1]:grid_spacing(newgrid)[1]:newgrid.size[2]-grid_spacing(newgrid)[1]
    newy = newgrid.size[1]:grid_spacing(newgrid)[2]:newgrid.size[2]-grid_spacing(newgrid)[2]
    field = set_field(newgrid, func)
    p = contour(newx, newy, field', levels=levels, color=color, clabels=clabels, cbar=cbar, lw=lw, legend=legend, aspect_ratio=aspect_ratio, kwargs...)
    if show_grid
        for i in eachindex(x)
            plot!(p, [x[i], x[i]], [y[1], y[end]], color=:black, lw=0.3)
            plot!(p, [x[1], x[end]], [y[i], y[i]], color=:black, lw=0.3)
        end
    end
    plot!(p, axis = axis, grid = plot_grid, xlims = (x[1], x[end]), ylims = (y[1], y[end]))
    return(p)
end

"""
Plot a 2D level-set countour
field_contour(grid, field)
"""
function field_contour(grid::Grid, field; show_grid=true, 
    levels=6,
    color=:turbo,
    clabels=false,
    cbar=false,
    lw=2,
    legend=false,
    aspect_ratio=:equal,
    axis = false,
    plot_grid = false,
    kwargs...)
    if length(grid.dims) != 2
        throw("The grid must be 2D.")
    end
    x = grid.size[1]:grid_spacing(grid)[1]:grid.size[2]-grid_spacing(grid)[1]
    y = grid.size[1]:grid_spacing(grid)[2]:grid.size[2]-grid_spacing(grid)[2]
    p = contour(x, y, field', levels=levels, color=color, clabels=clabels, cbar=cbar, lw=lw, legend=legend, aspect_ratio=aspect_ratio, kwargs...)
    if show_grid
        for i in eachindex(x)
            plot!(p, [x[i], x[i]], [y[1], y[end]], color=:black, lw=0.3)
            plot!(p, [x[1], x[end]], [y[i], y[i]], color=:black, lw=0.3)
        end
    end
    plot!(p, axis = axis, grid = plot_grid)
    return(p)
end

"""
Plot a 2D implicit curve, scatter points can be added.
implicit_plot_2D(grid, func)
"""
function implicit_plot_2D(grid::Grid, func; 
    show_grid = true,
    lw = 2,
    dataset = [0.0 0.0;0.0 0.0], 
    markersize = 0,
    axis = false,
    plot_grid = false,
    kwargs...
    )
    f(x,y) = func([x,y])
    if length(grid.dims) != 2
        throw("The grid must be 2D.")
    end
    x = grid.size[1]:grid_spacing(grid)[1]:grid.size[2]-grid_spacing(grid)[1]
    y = grid.size[1]:grid_spacing(grid)[2]:grid.size[2]-grid_spacing(grid)[2]

    p = implicit_plot(f; xlims=grid.size, ylims=grid.size, legend=false, lw=lw)
    if show_grid
        for i in eachindex(x)
            plot!(p, [x[i], x[i]], [y[1], y[end]], color=:black, lw=0.3)
            plot!(p, [x[1], x[end]], [y[i], y[i]], color=:black, lw=0.3)
        end
    end
    plot!(p, dataset[1,:], dataset[2,:], seriestype=:scatter, markersize=markersize, axis = axis, grid = plot_grid, kwargs...)
    return(p)
end