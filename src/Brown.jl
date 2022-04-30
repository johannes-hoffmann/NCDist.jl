###############################################################################
#
#   (Density of the) Brown measure via subordination
#   for arbitrary variable distributions at arbitrary points
#
###############################################################################

"""
    brown_rational(L, var_dists, eval_points; <keyword arguments>)

Compute the density of the Brown measure of a self-adjoint  rational
 expression `L` with respect to distributions of the variables `var_dists` at
 certain two-dimensional points `eval_points`.
Returns an array `z_values` such that `z_values[k]` is the value computed for
 the point `eval_points[k]`.

# Arguments

- `L::AbstractArray{Complex{Float64}, 3}`: the matrix-valued coefficients of a
  self-adjoint linearization of the rational expression.
- `var_dists::AbstractVector{<:Function}`: Cauchy transforms on the upper
  half-plane of the distributions of the variables.
- `eval_points::Vector{NTuple{2, <:Float64}}`: Real two-dimensional points to
  evaluate the (approximation of the density of) Brown measure at.
- `x_shift::Float64`, `y_shift::Float64`: step-sizes in x- resp. y-direction
  for computation of the derivative.

Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
 `L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function brown_rational(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    eval_points::Vector{NTuple{2, <:Float64}};
    kwargs...
)
    return brown_worker(L, var_dists, eval_points, false; kwargs...)
end

"""
    brown_operator(L, var_dists, eval_points; <keyword arguments>)

Compute the density of the Brown measure of an operator-valued element `L` with
 respect to distributions of the variables `var_dists` at certain
 two-dimensional points `eval_points`.
Returns an array `z_values` such that `z_values[k]` is the value computed for
 the point `eval_points[k]`.

# Arguments

- `L::AbstractArray{Complex{Float64}, 3}`: the matrix-valued coefficients of
  an operator-valued element.
- `var_dists::AbstractVector{<:Function}`: Cauchy transforms on the upper
  half-plane of the distributions of the variables.
- `eval_points::Vector{NTuple{2, <:Float64}}`: Real two-dimensional points to
  evaluate the (approximation of the density of) Brown measure at.
- `x_shift::Float64`, `y_shift::Float64`: step-sizes in x- resp. y-direction
  for computation of the derivative.

Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
 `L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function brown_operator(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    eval_points::Vector{NTuple{2, <:Float64}};
    kwargs...
)
    return brown_worker(L, var_dists, eval_points, true; kwargs...)
end

function brown_worker(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    eval_points::Vector{NTuple{2, <:Float64}},
    use_operator::Bool;
    x_shift::Float64 = 1.0e-3,
    y_shift::Float64 = 1.0e-3,
    kwargs...
)
    check_input(L, var_dists)
    if use_operator
        L = hermitize(L)
    end
    (eigen_val_tensor, eigen_vec_tensor, dim_range) = eigen_tensor(L)
    const_mat = L[:, :, 1]
    z_values = Float64[]
    if use_operator
        cauchy_brown = cauchy_brown_operator
    else
        cauchy_brown = cauchy_brown_rational
    end
    cauchy(z) = cauchy_brown(
        var_dists,
        size(L, 1), # matrix dimension
        size(L, 3), # number of variables + 1
        dim_range,
        const_mat,
        eigen_val_tensor,
        eigen_vec_tensor,
        z;
        kwargs...
    )
    for (x, y) in eval_points
        z_0 = cauchy(x + 1im * y)
        z_x = cauchy((x + x_shift) + 1im * y)
        z_y = cauchy(x + 1im * (y + y_shift))
        z_1 = (z_x - z_0) / x_shift
        z_2 = (z_y - z_0) / y_shift
        z = 1 / (2 * pi) * real(z_1 + 1im * z_2)
        push!(z_values, z)
    end
    return z_values
end
###############################################################################
#
#   (Density of the) Brown measure via subordination
#   for arbitrary variable distributions on a grid
#
###############################################################################

"""
    brown_rational_grid(L, var_dists, x_range, y_range)

Compute the density of the Brown measure of a self-adjoint rational expression
 `L` with respect to distributions of the variables `var_dists` on a grid
 specified by `x_range` and `y_range`.
Returns a matrix `z_values` such that `z_values[k, l]` is the value computed
 for the point `(x_range[k], y_range[l])`.

# Arguments

- `L::AbstractArray{Complex{Float64}, 3}`: the matrix-valued coefficients of a
  self-adjoint linearization of the rational expression.
- `var_dists::AbstractVector{<:Function}`: Cauchy transforms on the upper
  half-plane of the distributions of the variables.
- `x_range::AbstractRange{<:Float64}`, `y_range::AbstractRange{<:Float64}`: the
  Cartesian product of `x_range` and `y_range` is the set of all points to
   evaluate the (approximation of the density of) Brown measure at.

Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
 `L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function brown_rational_grid(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    x_range::AbstractRange{<:Float64},
    y_range::AbstractRange{<:Float64};
    kwargs...
)
    return brown_worker_grid(L, var_dists, x_range, y_range, false; kwargs...)
end

"""
    brown_operator_grid(L, var_dists, x_values, y_values)

Compute the density of the Brown measure of of an operator-valued element `L`
 with respect to distributions of the variables `var_dists` on a grid
 specified by `x_range` and `y_range`.
Returns a matrix `z_values` such that `z_values[k, l]` is the value computed
 for the point `(x_range[k], y_range[l])`.

# Arguments

- `L::AbstractArray{Complex{Float64}, 3}`: the matrix-valued coefficients of an
  operator-valued element.
- `var_dists::AbstractVector{<:Function}`: Cauchy transforms on the upper
  half-plane of the distributions of the variables.
- `x_range::AbstractRange{<:Float64}`, `y_range::AbstractRange{<:Float64}`: the
  Cartesian product of `x_range` and `y_range` is the set of all points to
   evaluate the (approximation of the density of) Brown measure at.

Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
 `L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function brown_operator_grid(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    x_range::AbstractRange{<:Float64},
    y_range::AbstractRange{<:Float64};
    kwargs...
)
    return brown_worker_grid(L, var_dists, x_range, y_range, true; kwargs...)
end

function brown_worker_grid(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    x_range::AbstractRange{<:Float64},
    y_range::AbstractRange{<:Float64},
    use_operator::Bool;
    kwargs...
)
    check_input(L, var_dists)
    if use_operator # Hermitization
        L = hermitize(L)
    end # end Hermitization
    (eigen_val_tensor, eigen_vec_tensor, dim_range) = eigen_tensor(L)
    const_mat = L[:, :, 1]
    Nx = length(x_range)
    Ny = length(y_range)
    dx = step(x_range)
    dy = step(y_range)
    G = zeros(Complex{Float64}, Nx + 1, Ny + 1)
    if use_operator
        cauchy = cauchy_brown_operator
    else
        cauchy = cauchy_brown_rational
    end
    for x_index in 1:(Nx + 1) # one step further than last(x_range)
        @debug "computed $x_index rows out of $(Nx + 1) for G in brown_operator_grid"
        for y_index in 1:(Ny + 1)  # one step further than last(y_range)
            x = first(x_range) + (x_index - 1) * dx
            y = first(y_range) + (y_index - 1) * dy
            G[x_index, y_index] = cauchy(
                var_dists,
                size(L, 1), # matrix dimension
                size(L, 3), # number of variables + 1
                dim_range,
                const_mat,
                eigen_val_tensor,
                eigen_vec_tensor,
                x + 1im * y;
                kwargs...
            )
        end
    end
    z_matrix = zeros(Float64, Nx, Ny)
    for x_index in 1:Nx, y_index in 1:Ny
        z_1 = (G[x_index + 1, y_index] - G[x_index, y_index]) / dx
        z_2 = (G[x_index, y_index + 1] - G[x_index, y_index]) / dy
        z_matrix[x_index, y_index] = 1 / (2 * pi) * real(z_1 + im * z_2)
    end
    return z_matrix
end

function cauchy_brown_operator(
    var_dists::Vector{<:Function},
    mat_dim::Int,
    nvars::Int,
    dim_range::Array{Int, 1},
    const_mat::Array{Complex{Float64}, 2},
    eigen_val_tensor::Array{Float64, 3},
    eigen_vec_tensor::Array{Complex{Float64}, 3},
    z::Complex{Float64};
    epsilon::Float64 = 1.0e-3,
    kwargs...
)
    b = Matrix(UniformScaling(1im * epsilon)(mat_dim))
    odim = div(mat_dim, 2)
    b[1:odim, (odim + 1):mat_dim] = UniformScaling(z)(odim)
    b[(odim + 1):mat_dim, 1:odim] = UniformScaling(conj(z))(odim)
    G = inv(F_transform_XY_brown(
        var_dists,
        mat_dim,
        nvars,
        dim_range,
        eigen_val_tensor,
        eigen_vec_tensor,
        nvars - 1,
        b - const_mat;
        kwargs...
    ))
    return tr(G[(odim + 1):mat_dim, 1:odim]) / odim
end

function cauchy_brown_rational(
    var_dists::Vector{<:Function},
    mat_dim::Int,
    nvars::Int,
    dim_range::Array{Int, 1},
    const_mat::Array{Complex{Float64}, 2},
    eigen_val_tensor::Array{Float64, 3},
    eigen_vec_tensor::Array{Complex{Float64}, 3},
    z::Complex{Float64};
    epsilon::Float64 = 1.0e-3,
    kwargs...
)
    b = Matrix(UniformScaling(1im * epsilon)(mat_dim))
    b[1, 2] = z
    b[2, 1] = conj(z)
    G = inv(F_transform_XY_brown(
        var_dists,
        mat_dim,
        nvars,
        dim_range,
        eigen_val_tensor,
        eigen_vec_tensor,
        nvars - 1,
        b - const_mat;
        kwargs...
    ))
    return G[2, 1]
end

###############################################################################
#
#   (Density of the) Brown measure with semicircle variable distributions
#   at arbitrary points
#
###############################################################################

"""
    brown_rational_semicircle(L, eval_points; <keyword arguments>)

Compute the density of the Brown measure of a self-adjoint rational
 expression `L` with respect to semicircle distributions of the variables at
 certain two-dimensional points `eval_points`.
Returns an array `z_values` such that `z_values[k]` is the value computed for
 the point `eval_points[k]`.

# Arguments

- `L::AbstractArray{Complex{Float64}, 3}`: the matrix-valued coefficients of a
  self-adjoint linearization of the rational expression.
- `eval_points::Vector{NTuple{2, <:Float64}}`: Real two-dimensional points to
  evaluate the (approximation of the density of) Brown measure at.
- `x_shift::Float64`, `y_shift::Float64`: step-sizes in x- resp. y-direction
  for computation of the derivative.

Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
 `L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function brown_rational_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    eval_points::Vector{NTuple{2, <:Float64}};
    kwargs...
)
    return brown_worker_semicircle(L, eval_points, false; kwargs...)
end

"""
    brown_operator_semicircle(L, eval_points; <keyword arguments>)

Compute the density of the Brown measure of an operator-valued element `L` with
 respect to semicircle distributions of the variables at certain
 two-dimensional points `eval_points`.
Returns an array `z_values` such that `z_values[k]` is the value computed for
 the point `eval_points[k]`.

# Arguments

- `L::AbstractArray{Complex{Float64}, 3}`: the matrix-valued coefficients of
  an operator-valued element.
- `eval_points::Vector{NTuple{2, <:Float64}}`: Real two-dimensional points to
  evaluate the (approximation of the density of) Brown measure at.
- `x_shift::Float64`, `y_shift::Float64`: step-sizes in x- resp. y-direction
  for computation of the derivative.

Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
 `L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function brown_operator_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    eval_points::Vector{NTuple{2, <:Float64}};
    kwargs...
)
    return brown_worker_semicircle(L, eval_points, true; kwargs...)
end

function brown_worker_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    eval_points::Vector{NTuple{2, <:Float64}},
    use_operator::Bool;
    x_shift::Float64 = 1.0e-3,
    y_shift::Float64 = 1.0e-3,
    epsilon::Float64 = 1.0e-3,
    kwargs...
)
    check_input(L)
    if use_operator
        L = hermitize(L)
    end
    z_values = Float64[]
    cauchy(z) = cauchy_brown_semicircle(L, z, epsilon, use_operator; kwargs...)
    for (x, y) in eval_points
        z_0 = cauchy(x + 1im * y)
        z_x = cauchy((x + x_shift) + 1im * y)
        z_y = cauchy(x + 1im * (y + y_shift))
        z_1 = (z_x - z_0) / x_shift
        z_2 = (z_y - z_0) / y_shift
        z = 1 / (2 * pi) * real(z_1 + 1im * z_2)
        push!(z_values, z)
    end
    return z_values
end

###############################################################################
#
#   (Density of the) Brown measure with semicircle variable distributions
#   on a grid
#
###############################################################################

"""
    brown_rational_grid_semicircle(L, x_values, y_values)

Compute the density of the Brown measure of a self-adjoint rational expression
 `L` with respect to semicircle distributions of the variables on a grid
 specified by `x_range` and `y_range`.
Returns a matrix `z_values` such that `z_values[k, l]` is the value computed
 for the point `(x_range[k], y_range[l])`.

# Arguments

- `L::AbstractArray{Complex{Float64}, 3}`: the matrix-valued coefficients of a
  self-adjoint linearization of the rational expression.
- `x_range::AbstractRange{<:Float64}`, `y_range::AbstractRange{<:Float64}`: the
  Cartesian product of `x_range` and `y_range` is the set of all points to
   evaluate the (approximation of the density of) Brown measure at.

Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
 `L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function brown_rational_grid_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    x_range::AbstractRange{<:Float64},
    y_range::AbstractRange{<:Float64};
    kwargs...
)
    return brown_worker_grid_semicircle(L, x_range, y_range, false; kwargs...)
end

"""
    brown_operator_grid_semicircle(L, x_values, y_values)

Compute the density of the Brown measure of of an operator-valued element `L`
 with respect to to semicircle distributions of the variables on a grid
 specified by `x_range` and `y_range`.
Returns a matrix `z_values` such that `z_values[k, l]` is the value computed
 for the point `(x_range[k], y_range[l])`.

# Arguments

- `L::AbstractArray{Complex{Float64}, 3}`: the matrix-valued coefficients of an
  operator-valued element.
- `x_range::AbstractRange{<:Float64}`, `y_range::AbstractRange{<:Float64}`: the
  Cartesian product of `x_range` and `y_range` is the set of all points to
   evaluate the (approximation of the density of) Brown measure at.

Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
 `L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function brown_operator_grid_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    x_range::AbstractRange{<:Float64},
    y_range::AbstractRange{<:Float64};
    kwargs...
)
    return brown_worker_grid_semicircle(L, x_range, y_range, true; kwargs...)
end

function brown_worker_grid_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    x_range::AbstractRange{<:Float64},
    y_range::AbstractRange{<:Float64},
    use_operator::Bool;
    epsilon::Float64 = 1.0e-3,
    kwargs...
)
    check_input(L)
    if use_operator
        L = hermitize(L)
    end
    Nx = length(x_range)
    Ny = length(y_range)
    dx = step(x_range)
    dy = step(y_range)
    G = zeros(Complex{Float64}, Nx + 1, Ny + 1)
    cauchy(z) = cauchy_brown_semicircle(L, z, epsilon, use_operator; kwargs...)
    for x_index in 1:(Nx + 1) # one step further than last(x_range)
        @debug "computed $x_index rows out of $(Nx + 1) for G in brown_operator_grid"
        for y_index in 1:(Ny + 1)  # one step further than last(y_range)
            x = first(x_range) + (x_index - 1) * dx
            y = first(y_range) + (y_index - 1) * dy
            G[x_index, y_index] = cauchy(x + 1im * y)
        end
    end
    z_matrix = zeros(Float64, Nx, Ny)
    for x_index in 1:Nx, y_index in 1:Ny
        z_1 = (G[x_index + 1, y_index] - G[x_index, y_index]) / dx
        z_2 = (G[x_index, y_index + 1] - G[x_index, y_index]) / dy
        z_matrix[x_index, y_index] = 1 / (2 * pi) * real(z_1 + im * z_2)
    end
    return z_matrix
end

function cauchy_brown_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    z::Complex{Float64},
    epsilon::Float64,
    use_operator::Bool;
    kwargs...
)
    mat_dim = size(L, 1)
    if use_operator
        b = Matrix(UniformScaling(1im * epsilon)(mat_dim))
        odim = div(mat_dim, 2)
        b[1:odim, (odim + 1):mat_dim] = UniformScaling(z)(odim)
        b[(odim + 1):mat_dim, 1:odim] = UniformScaling(conj(z))(odim)
    else
        b = Matrix(UniformScaling(1im * epsilon)(mat_dim))
        b[1, 2] = z
        b[2, 1] = conj(z)
    end
    W = cauchy_operator_semicircle_iteration(L, b; kwargs...)
    if use_operator
        return tr(W[(odim + 1):mat_dim, 1:odim]) / odim
    else
        return W[2, 1]
    end
end
