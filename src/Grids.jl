"""
structure to organize grid and interface information
currently we only consider uniform cubic grid
"""
struct Grid
    dims
    size::Tuple{Float64, Float64}
    len::Float64
end

function Grid(dims, size)
    len = abs(size[2] - size[1])
    return Grid(dims, size, len)
end

"""
     grid_spacing(grid::Grid)

Return a Vector containing the grid spacing in each dimension
"""
function grid_spacing(grid::Grid)
    if length(grid.dims) == 2
        return [grid.len/grid.dims[1], grid.len/grid.dims[2]]
    elseif  length(grid.dims) == 3
        return [grid.len/grid.dims[1], grid.len/grid.dims[2], grid.len/grid.dims[3]]
    else
        throw("Grid dimension not supported")
    end
end

"""
    pos(grid::Grid, i)

Return the position of an index
"""
function pos(grid::Grid, i)
    if length(grid.dims) == 2
        return Float64[grid.size[1] + grid.len*(i[1] - 1)/grid.dims[1], grid.size[1] + grid.len*(i[2] - 1)/grid.dims[2]]
    elseif  length(grid.dims) == 3
        return Float64[grid.size[1] + grid.len*(i[1] - 1)/grid.dims[1], grid.size[1] + grid.len*(i[2] - 1)/grid.dims[2], grid.size[1] + grid.len*(i[3] - 1)/grid.dims[3]]
    else
        throw("Grid dimension not supported")
    end
end

"""
    set_field(grid::Grid, func::Function)

Get a 2D field for the `grid` from a function f.
"""
function set_field(grid::Grid, func::Function)
    field = CircularArray(0.0, grid.dims)
    for I in eachindex(field)
        x = pos(grid, Tuple(I))
        field[I] = func(x)
    end
    return field
end

"""
Create a vector representing the grid field
"""
function field_to_vector(field, grid::Grid)
    sz = grid.dims
    N  = sz[1]*sz[2]
    out = zeros(N)
    n = 1
    for I in eachindex(field)
        out[n] = field[I]
        n = n + 1
    end
    return out
end

"""
Create a field representing the vector
"""
function vector_to_field(vector, grid::Grid)
    field = CircularArray(0.0, grid.dims)
    n = 1
    for I in eachindex(field)
        field[I] = vector[n]
        n = n + 1
    end
    return field
end

"""
find index of a point [i,j] in the vector
"""
function find_index(label, grid::Grid)
    sz = grid.dims
    return (label[2] - 1)*sz[2] + label[1]
end


