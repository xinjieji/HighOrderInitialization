# This file includes functions that is useful for calculating a evolving surface
# The basic function for evolving surface: dϕ/dt = -v·∇ϕ

# shorthands for adding entries to a sparse matrix
function add_entry(i::Int64, j::Int64, v::Float64, I::Vector{Int64}, J::Vector{Int64}, V::Vector{Float64})
    push!(I, i)
    push!(J, j)
    push!(V, v)
end

"""
Given the grid, generate the first derivative operators.
"""
function create_first_derivative_matrix(grid::Grid)
    h = grid_spacing(grid)
    if h[1] != h[2]
        throw("Grid spacing not equal in x and y directions")
    end

    # centered
    ind0 = -1:1#-2:2 #-1:1
    s0 = OffsetArray([-1/2, 0, 1/2]/h[1], ind0)#OffsetArray([1/12, -8/12, 0., 8/12, -1/12]/h[1], ind0)#OffsetArray([-1/2, 0, 1/2]/h[1], ind0)
    # forward
    ind1 = 0:2
    s1 = OffsetArray([-3/2, 2, -1/2]/h[1], ind1)
    # backward
    ind_1 = -2:0
    s_1 = OffsetArray([1/2, -2, 3/2]/h[1], ind_1)

    Ng = grid.dims[1]*grid.dims[2]

    # x-direction
    I_x, J_x, V_x = Int64[], Int64[], Float64[]
    row  = 1
    for j in 1:grid.dims[2]
        for i in 1:grid.dims[1]
            if i == 1 #|| i == 2
                for di in ind1
                    col = find_index([i+di,j],grid)
                    add_entry(row, col, s1[di], I_x, J_x, V_x)
                end
            elseif i == grid.dims[1] #|| i == grid.dims[1]-1
                for di in ind_1
                    col = find_index([i+di,j],grid)
                    add_entry(row, col, s_1[di], I_x, J_x, V_x)
                end
            else
                for di in ind0
                    col = find_index([i+di,j],grid)
                    add_entry(row, col, s0[di], I_x, J_x, V_x)
                end
            end
            row += 1
        end
    end
    Dx = sparse(I_x, J_x, V_x, Ng, Ng)

    # y-direction
    I_y, J_y, V_y = Int64[], Int64[], Float64[]
    row  = 1
    for j in 1:grid.dims[1]
        for i in 1:grid.dims[2]
            if j == 1 #|| j == 2
                for dj in ind1
                    col = find_index([i,j+dj],grid)
                    add_entry(row, col, s1[dj], I_y, J_y, V_y)
                end
            elseif j == grid.dims[2] #|| j == grid.dims[2]-1
                for dj in ind_1
                    col = find_index([i,j+dj],grid)
                    add_entry(row, col, s_1[dj], I_y, J_y, V_y)
                end
            else
                for dj in ind0
                    col = find_index([i,j+dj],grid)
                    add_entry(row, col, s0[dj], I_y, J_y, V_y)
                end
            end
            row += 1
        end
    end
    Dy = sparse(I_y, J_y, V_y, Ng, Ng)

    return Dx, Dy
end