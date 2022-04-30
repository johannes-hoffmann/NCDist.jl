###############################################################################
#
#   Distributions via subordination with arbitrary variable distributions
#
###############################################################################

"""
    distribution_rational(L, var_dists, x_values)

Compute the distribution of a (self-adjoint) rational expression, given via the
matrix-valued coefficients `L` of a self-adjoint linearization, with respect to
 distributions of the variables, given via their Cauchy transforms (on the
 upper half-plane) `var_dists`, at certain real points, given via `x_values`.
Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
`L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function distribution_rational(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    x_values::AbstractVector{<:Float64};
    kwargs...
)
    return distribution_worker(L, var_dists, x_values, false; kwargs...)
end

"""
    distribution_operator(L, var_dists, x_values)

Compute the distribution of an operator-valued element, given via the
self-adjoint matrix-valued coefficients `L`, with respect to distributions of
the variables, given via their Cauchy transforms (on the upper half-plane)
`var_dists`, at certain real points, given via `x_values`.
Note that `L[:, :, 1]` corresponds to the constant coefficients, thus
`L[:, :, n]` corresponds to the variable with distribution `var_dists[n - 1]`.
"""
function distribution_operator(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    x_values::AbstractVector{<:Float64};
    kwargs...
)
    return distribution_worker(L, var_dists, x_values, true; kwargs...)
end

function distribution_worker(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function},
    x_values::AbstractVector{<:Float64},
    use_operator::Bool;
    y_shift::Float64 = 1.0e-6,
    kwargs...
)
    check_input(L, var_dists)
    mat_dim = size(L, 1)
    nvars = size(L, 3) - 1
    if y_shift <= 0
        error("y_shift must be strictly positive, but is $y_shift")
    end
    (eigen_val_tensor, eigen_vec_tensor, dim_range) = eigen_tensor(L)
    const_mat = L[:, :, 1]
    y_values = Float64[]
    for x_value in x_values
        z = x_value + y_shift * im
        if use_operator
            b = UniformScaling(z)(mat_dim)
        else
            b = UniformScaling(y_shift * im)(mat_dim)
            b[1, 1] = z
        end
        G = inv(F_transform_XY_dist(
            var_dists,
            mat_dim,
            dim_range,
            eigen_val_tensor,
            eigen_vec_tensor,
            nvars,
            b - const_mat;
            approximation_error = y_shift/10,
            kwargs...
        ))
        if use_operator
            temp_value = tr(G) / mat_dim
        else
            temp_value = G[1, 1]
        end
        push!(y_values, -imag(temp_value) / pi)
    end
    return y_values
end

###############################################################################
#
#   Distributions with semicircle variable distributions
#
###############################################################################

"""
    distribution_rational_semicircle(L, var_dists, x_values)

Compute the distribution of a (self-adjoint) rational expression, given via the
matrix-valued coefficients `L` of a self-adjoint linearization, with respect to
 semicircle distributions of the variables at certain real points, given via
 `x_values`.
Note that `L[:, :, 1]` corresponds to the constant coefficients.
"""
function distribution_rational_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    x_values::AbstractVector{<:Float64};
    kwargs...
)
    return distribution_worker_semicircle(L, x_values, false; kwargs...)
end

"""
    distribution_operator_semicircle(L, var_dists, x_values)

Compute the distribution of an operator-valued element, given via the
self-adjoint matrix-valued coefficients `L`, with respect to semicircle
distributions of the variables at certain real points, given via `x_values`.
Note that `L[:, :, 1]` corresponds to the constant coefficients.
"""
function distribution_operator_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    x_values::AbstractVector{<:Float64};
    kwargs...
)
    return distribution_worker_semicircle(L, x_values, true; kwargs...)
end

function distribution_worker_semicircle(
    L::AbstractArray{Complex{Float64}, 3},
    x_values::AbstractVector{<:Float64},
    use_operator::Bool;
    y_shift::Float64 = 1.0e-6,
    kwargs...
)
    check_input(L)
    mat_dim = size(L, 1)
    y_values = Float64[]
    for x_value in x_values
        if use_operator
            b = UniformScaling(x_value + y_shift * im)(mat_dim)
        else
            b = UniformScaling(y_shift * im)(mat_dim)
            b[1, 1] = x_value + y_shift * im
        end
        G = cauchy_operator_semicircle_iteration(L, b; kwargs...)
        if use_operator
            temp_value = tr(G) / mat_dim
        else
            temp_value = G[1, 1]
        end
        push!(y_values, -imag(temp_value) / pi)
    end
    return y_values
end
