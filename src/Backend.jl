###############################################################################
#
#   Subordination backend
#
###############################################################################

function F_transform_XY_dist(
    var_dists::AbstractVector{<:Function},
    mat_dim::Int,
    dim_range::Vector{Int},
    eigen_val_tensor::Array{Float64, 3},
    eigen_vec_tensor::Array{Complex{Float64}, 3},
    step::Int,
    b::AbstractMatrix{Complex{Float64}};
    iteration_maximum::Int = 10_000,
    approximation_error::Float64 = 1.0e-3
)
    if step == 1
        return F_transform_X(var_dists, 1, mat_dim, dim_range[1],
            eigen_val_tensor[:, :, 1], eigen_vec_tensor[:, :, 1], b
        )
    end
    sub = UniformScaling(1.0im)(mat_dim)
    F = F_transform_X(var_dists, step, mat_dim, dim_range[step],
        eigen_val_tensor[:, :, step], eigen_vec_tensor[:, :, step], sub
    )
    new_sub = F_transform_XY_dist(var_dists, mat_dim, dim_range,
        eigen_val_tensor, eigen_vec_tensor, step - 1, b + F - sub
    ) - F + sub
    exitcounter = 1
    while exitcounter < iteration_maximum && norm(sub - new_sub) > approximation_error # Frobenius norm
        exitcounter += 1
        sub = new_sub
        F = F_transform_X(var_dists, step, mat_dim, dim_range[step],
            eigen_val_tensor[:, :, step], eigen_vec_tensor[:, :, step], sub
        )
        new_sub = F_transform_XY_dist(var_dists, mat_dim, dim_range,
            eigen_val_tensor, eigen_vec_tensor, step - 1, b + F - sub
        ) - F + sub
    end
    if exitcounter == iteration_maximum
        @warn "too many iterations needed"
    end
    return F_transform_X(var_dists, step, mat_dim, dim_range[step],
        eigen_val_tensor[:, :, step], eigen_vec_tensor[:, :, step], sub
    )
end

function F_transform_XY_brown(
    var_brown::AbstractVector{<:Function},
    dim::Int,
    nvars::Int,
    dim_range::Vector{Int},
    eigen_val_tensor::Array{Float64, 3},
    eigen_vec_tensor::Array{Complex{Float64}, 3},
    step::Int,
    b::AbstractMatrix{Complex{Float64}};
    approximation_error::Float64 = 1.0e-3
)
    save = 50
    average = 10
    exit = 10
    if step == 1
        return F_transform_X(var_brown, 1, dim, dim_range[1],
            eigen_val_tensor[:, :, 1], eigen_vec_tensor[:, :, 1], b
        )
    end
    sub = UniformScaling(1.0im)(dim)
    F = F_transform_X(var_brown, step, dim, dim_range[step],
        eigen_val_tensor[:, :, step], eigen_vec_tensor[:, :, step], sub
    )
    new_sub = F_transform_XY_brown(var_brown, dim, nvars, dim_range,
        eigen_val_tensor, eigen_vec_tensor, step - 1, b + F - sub
    ) - F + sub
    save_sub = zeros(Complex{Float64}, size(new_sub, 1), size(new_sub, 2), save)
    save_sub[:, :, 1] = new_sub
    exitcounter = 1
    counter = 1
    while exitcounter < exit && norm(sub - new_sub) > approximation_error # Frobenius norm
        if counter < save
            counter = counter + 1
            sub = new_sub
            F = F_transform_X(var_brown, step, dim, dim_range[step],
                eigen_val_tensor[:, :, step], eigen_vec_tensor[:, :, step], sub
            )
            new_sub = F_transform_XY_brown(var_brown, dim, nvars, dim_range,
                eigen_val_tensor, eigen_vec_tensor, step - 1, b + F - sub
            ) - F + sub
            save_sub[:, :, counter] = new_sub
        else
            counter = 1
            exitcounter = exitcounter + 1
            new_sub = zeros(dim, dim)
            for k in 1:average
                new_sub = new_sub + save_sub[:, :, save + 1 - k]
            end
            new_sub = (1 / average) * new_sub
            save_sub[:, :, 1] = new_sub
        end
    end
    if exitcounter == exit
        @warn "too many iterations needed"
    end
    return F_transform_X(var_brown, step, dim, dim_range[step],
        eigen_val_tensor[:, :, step], eigen_vec_tensor[:, :, step], new_sub
    )
end

function F_transform_X(
    var_dists::AbstractVector{<:Function},
    j::Int,
    mat_dim::Int,
    d::Int,
    eigen_val_mat::Matrix{Float64},
    eigen_vec_mat::Matrix{Complex{Float64}},
    b::AbstractMatrix{Complex{Float64}}
)
    beta = adjoint(eigen_vec_mat) * b * eigen_vec_mat
    N = size(beta, 1)
    A = beta[1:d, 1:d]
    B = beta[1:d, (d + 1):N]
    C = beta[(d + 1):N, 1:d]
    D = beta[(d + 1):N, (d + 1):N]
    # A B
    # C D
    mu = eigen_val_mat[1:d, 1:d]
    S = A - B * inv(D) * C
    Eig = eigen(inv(mu) * S)
    F = zeros(Complex{Float64}, d, d)
    for k in 1:d
        z = complex(Eig.values[k])
        f = var_dists[j]
        if imag(z) < 0
            f_z = conj(f(conj(z)))
        elseif imag(z) == 0
            f_z = f(z)
        else
            f_z = f(z)
        end
        F[k, k] = 1 / f_z
    end
    F = mu * Eig.vectors * F * inv(Eig.vectors)
    Q1 = eigen_vec_mat * Schur2(B * inv(D))
    Q2 = Schur1(inv(D) * C) * adjoint(eigen_vec_mat)
    y = [F zeros(d, mat_dim - d); zeros(mat_dim - d, d) D]
    return Q1 * y * Q2
end

function Schur1(beta::AbstractMatrix{Complex{Float64}})
    (p, q) = size(beta)
    return [
        Matrix{Complex{Float64}}(I, q, q) zeros(q, p);
        beta Matrix{Complex{Float64}}(I, p, p)
    ]
end

function Schur2(beta::AbstractMatrix{Complex{Float64}})
    (p, q) = size(beta)
    return [
        Matrix{Complex{Float64}}(I, p, p) beta;
        zeros(q, p) Matrix{Complex{Float64}}(I, q, q)
    ]
end

function eigen_tensor(L::AbstractArray{Complex{Float64}, 3})
    mat_dim = size(L, 1)
    nvars = size(L, 3) - 1
    eigen_val_tensor = zeros(Float64, mat_dim, mat_dim, nvars)
    eigen_vec_tensor = zeros(Complex{Float64}, mat_dim, mat_dim, nvars)
    dim_range = zeros(Int, nvars)
    for var_index in 1:nvars
        (vals, vecs) = sorted_eigen(Hermitian(L[:, :, var_index + 1]))
        eigen_val_tensor[:, :, var_index] = Diagonal(real.(vals))
        eigen_vec_tensor[:, :, var_index] = vecs
        d = 0
        for l in 1:mat_dim
            if abs(eigen_val_tensor[l, l, var_index]) >= 0.000001
                d += 1
            end
        end
        dim_range[var_index] = d
    end
    return (eigen_val_tensor, eigen_vec_tensor, dim_range)
end

function sorted_eigen(A::AbstractMatrix{<:Number})
    n = size(A, 1)
    if n != size(A, 2)
        error("matrix is not quadratic")
    end
    (eig_vals, eig_vecs) = eigen(A)
    for k in 1:n
        for l in (k + 1):n
            if abs(eig_vals[l]) > abs(eig_vals[k])
                eig_vals[l], eig_vals[k] = eig_vals[k], eig_vals[l]
                eig_vecs[:, l], eig_vecs[:, k] = eig_vecs[:, k], eig_vecs[:, l]
            end
        end
    end
    return (eig_vals, eig_vecs)
end

###############################################################################
#
#   Semicircle iteration backend
#
###############################################################################

function cauchy_operator_semicircle_iteration(
    L::AbstractArray{Complex{Float64}, 3},
    b::AbstractMatrix{Complex{Float64}};
    iteration_const::Union{Int, Float64} = 0.5,
    iteration_maximum::Int = 200_000,
    additional_iterations::Int = 0,
    approximation_error::Float64 = 0.1, # error in θ # should be about 1/(4n)
    loop_breaking::Bool = true
)
    if iteration_const > 1 || iteration_const <= 0
        error("iteration_const has to be in the interval (0,1]")
    end
    imag_b = (b - adjoint(b)) / (2im)
    y_shift = inv(opnorm(inv(imag_b)))
    mat_dim = size(b, 1)
    approximation_error = minimum([approximation_error, 1 / (4 * mat_dim)])
    W = UniformScaling(-1.0im)(mat_dim)
    exitcounter = 0 # debug
    σ = y_shift * approximation_error / (1 + y_shift * approximation_error)
    # delta_bound = y_shift * approximation_error / (y_shift + approximation_error)
    delta_bound = σ * y_shift
    for iteration in 1:iteration_maximum
        if delta(L, b, W) <= delta_bound
            break
        end
        if loop_breaking && mod(iteration, 100) == 0 # random variation to break deadly loops
            random_const = rand(0.01:0.01:1)
            W = (1 - random_const) * W + random_const * inv(b - eta(L, W))
        else
            W = (1 - iteration_const) * W + iteration_const * inv(b - eta(L, W))
        end
        exitcounter += 1 # debug
    end
    for iteration in 1:additional_iterations
        W = (1 - iteration_const) * W + iteration_const * inv(b - eta(L, W))
    end
    d_str = ""
    if iteration_const != 0.5
        d_str *= "($(1 - iteration_const), $iteration_const): "
    end
    pad_len = length(string(iteration_maximum))
    d_str *= "$(lpad(string(exitcounter), pad_len, ' ')) of $iteration_maximum iterations needed, "
    d_str *= "delta: $(delta(L, b, W)) of $delta_bound, y_shift: $y_shift"
    if exitcounter == iteration_maximum
        @warn d_str
    else
        @debug d_str
    end
    return W # G(b)
    # b * G(b) = 1 + η(G(b)) * G(b)
end

function delta(
    L::AbstractArray{Complex{Float64}, 3},
    b::AbstractMatrix{Complex{Float64}},
    w_0::AbstractMatrix{Complex{Float64}}
)
    return opnorm(b - inv(w_0) - eta(L, w_0))
end

function eta(
    L::AbstractArray{Complex{Float64}, 3},
    b::AbstractMatrix{Complex{Float64}}
)
    n = size(L, 1)
    temp = zeros(Complex{Float64}, n, n)
    for i in 2:size(L, 3)
        temp += L[:, :, i] * b * adjoint(L[:, :, i])
    end
    return temp
end

###############################################################################
#
#   General backend
#
###############################################################################

# function continuation(f::Function)
#     return z -> if imag(z) < 0
#         conj(f(conj(z)))
#     else
#         f(z)
#     end
# end

function hermitize(L::AbstractArray{Complex{Float64}, 3})
    mat_dim = size(L, 1)
    nvars = size(L, 3) - 1
    ndim = 2 * mat_dim
    LHerm = zeros(Complex{Float64}, ndim, ndim, nvars + 1)
    for j in 1:(nvars + 1)
        LHerm[1:mat_dim, (mat_dim + 1):ndim, j] = L[:, :, j]
        LHerm[(mat_dim + 1):ndim, 1:mat_dim, j] = adjoint(L[:, :, j])
    end
    return LHerm
end

function check_input(L::AbstractArray{Complex{Float64}, 3})
    if size(L, 1) != size(L, 2)
        error("linearization matrices are not quadratic: $(size(L, 1)) rows
            vs $(size(L, 2)) columns")
    end
end

function check_input(
    L::AbstractArray{Complex{Float64}, 3},
    var_dists::AbstractVector{<:Function}
)
    check_input(L)
    if length(var_dists) != size(L, 3) - 1
        error("the number of variable distributions ($(length(var_dists))) does
            not match the number of variables ($(size(L, 3) - 1))")
    end
    return
end
