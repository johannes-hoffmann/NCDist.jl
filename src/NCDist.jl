module NCDist

import LinearAlgebra: Diagonal, eigen, Hermitian, I, UniformScaling
import LinearAlgebra: norm, opnorm, tr

export distribution_rational, distribution_operator
export distribution_rational_semicircle, distribution_operator_semicircle
export brown_rational, brown_rational_grid
export brown_operator, brown_operator_grid
export brown_rational_semicircle, brown_rational_grid_semicircle
export brown_operator_semicircle, brown_operator_grid_semicircle

export semicircle, poisson, bernoulli, symmetric_bernoulli, cauchy, arcsine

include("Backend.jl")
include("CauchyTransforms.jl")
include("Distribution.jl")
include("Brown.jl")

end # module