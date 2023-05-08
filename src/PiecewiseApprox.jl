# Currently the piecewise approximation uses the Taylor degree polynomials.

"""
This is the dictionary of the polynomial type.
From Saye(2014)

d -> dimension
nc -> number of coefficients
st -> number of points in the stencil
p -> expected order of accuracy

Dict((d,p) = (nc,st)
"""
const poly_dict_ = Dict(
    (2,3) => (6,12),
    (2,4) => (10,12),
    (2,5) => (15,24),
    (2,6) => (21,24),
    (3,3) => (10,32),
    (3,4) => (20,32),
    (3,5) => (35,88),
    (3,6) => (56,88)
)

"""
Create the Vandermonde matrix V for the given points.
x = [x_1 x_2 ... x_st], points with their positions.

For example, for 2D 3rd order polynomial, every line of V:
x_i -> [1 y y^2 x xy x^2]
"""
function vander(x, p::Int64)
    (d,np) = size(x) # d: dimension, np: number of points
    (nc,st) = poly_dict_[(d,p)] 

    v = zeros(np,nc)
    if d==2
        i = 1
        for dg1 = 0:p-1
            for dg2 = 0:p-1-dg1
                v[:,i] = (x[1,:].^dg1 .* x[2,:].^dg2)'
                i += 1
            end
        end
    elseif d==3 # Todo: need further test
        i = 1
        for dg1 = 0:p-1
            for dg2 = 0:p-1-dg1
                for dg3 = 0:p-1-dg1-dg2
                    v[:,i] = (x[1,:].^dg1 .* x[2,:].^dg2 .* x[3,:].^dg3)'
                    i += 1
                end
            end
        end
    else
        throw("Dimension not supported")
    end

    return v
end

"""
1. Detect cells that contain the interface.
2. Collect points for polynomials of the cells.
If the signs of ϕ at all vertices of the cell are not the same, 
then it covers the interface.
Note: at most we need two layers of points outside the cell

cell_index, X, ϕ_Vector = collect_points_2D(grid,ϕ,order)
- ϕ: the field of level-set
- cell_index: collection of the cells which cover the interface, labeled by the lower-left point index
- X: collection of the positions of points for the polynomials
- ϕ_Vector: collection of the values of ϕ at the points

Note: all the collections are based on the cell_index order, so X and ϕ_Vector will have duplicated entries.
"""
function collect_points_2D(grid::Grid, ϕ, p::Int64)
    d = length(grid.dims)
    if d != 2
        throw("Dimension not supported")
    end
    (nc,st) = poly_dict_[(d,p)]

    # 1. loop over all cells, find those covering the interface
    cell_index = [] # use the lower-left point index 
    for I in eachindex(ϕ)
        if I[1] == grid.dims[1] || I[2] == grid.dims[2]
            continue
        end
        if ϕ[I]*ϕ[I[1]+1, I[2]] < 0 || ϕ[I]*ϕ[I[1], I[2]+1] < 0 || ϕ[I]*ϕ[I[1]+1, I[2]+1] < 0
            push!(cell_index, I)
        end
    end

    # 2. collect points
    X = []
    ϕ_Vector = []
    for I in cell_index
        list = []
        if st==12
            # index distance list
            list = [(-1,1), (-1,0), (0,2), (0,1), (0,0), (0,-1),
            (1,2), (1,2), (1,0), (1,-1), (2,1), (2,0)]
        elseif st==24
            # index distance list
            list = [(i,j) for i in -1:2 for j in -1:2]
            push!(list, (-2,0), (-2,1), (3,0), (3,1),
                (0,-2), (0,3), (1,-2), (1,3))
        else
            throw("Stencil not supported")
        end
        # write the positions
        x = zeros(d, st)
        ϕ_vector = zeros(st)
        for i in 1:lastindex(list)
            x[:,i] = pos(grid, (I[1]+list[i][1], I[2]+list[i][2]))
            ϕ_vector[i] = ϕ[I[1]+list[i][1], I[2]+list[i][2]]
        end
        push!(X, x)
        push!(ϕ_Vector, ϕ_vector)
    end

    # 3. collect all the Vandermonde matrices
    V = []
    for x in X
        v = vander(x, p)
        if rank(v) != nc
            throw("The points are not linear independent")
        end
        push!(V, v)
    end

    return cell_index, X, ϕ_Vector, V
end

"""
Wrap collect_points_2D for scatter plot

Example:
> cell_index, X, ϕ_Vector = collect_points_2D(grid,ϕ,3)
> data = x_wrapper(X)
> implicit_plot_2D(grid,f, dataset=data, markersize = 1)
"""
function x_wrapper(X)
    data = X[1]
    for i = 2:lastindex(X)
        data = hcat(data, X[i])
    end
    return data
end

"""
Calculate the coefficients of the polynomial.
C is a matrix that stores all the coefficients, every column is a cell.
"""
function cal_coeff_2D(grid::Grid, ϕ, p::Int64)
    cell_index, X, ϕ_Vector, V = collect_points_2D(grid, ϕ, p)
    (st,nc) = size(V[1])
    C = zeros(nc, length(cell_index))
    for i in 1:lastindex(cell_index)
        C[:,i] = V[i] \ ϕ_Vector[i]
    end
    return cell_index, X, C
end

"""
Function of the interpolated polynomials in the cell, given the cell index and a point inside the cell,
output the approximated ϕ.
"""
function f_inter_2D(x, cell_num, C, p::Int64)
    c = C[:,cell_num]
    ϕ = 0
    i = 1
    for dg1 = 0:p-1
        for dg2 = 0:p-1-dg1
            ϕ += c[i]*(x[1]^dg1 * x[2]^dg2)
            i += 1
        end
    end
    return ϕ
end
