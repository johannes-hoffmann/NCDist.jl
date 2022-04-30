###############################################################################
#
#   Library of Cauchy transforms for important distributions
#
###############################################################################

###############################################################################
#
#   Semicircle distribution
#
###############################################################################

function semicircle(a::Real, t::Real)
    return semicircle(; a, t)
end

"""
    semicircle(a, t)
    semicircle(; a, t)

Return a function mapping a complex value `z` to the Cauchy transform of a
semicircle distribution with mean `a` and variance `t>0` at `z`
(default `a=0`, `t=1`).
"""
function semicircle(; a::Real = 0.0, t::Real = 1.0)
    if t <= 0
        error("variance t has to be positive, but is $t")
    end
    function f(z::Number)
        z = complex(z)
        return (z - a)/(2 * t) * (1 - sqrt(1 - 4 * t / (z - a)^2))
    end
    return f
end

###############################################################################
#
#   free Poisson distribution
#
###############################################################################

function poisson(lambda::Real, alpha::Real)
    return poisson(; lambda, alpha)
end

"""
    poisson(lambda, alpha)
    poisson(; lambda, alpha)

Return a function mapping a complex value `z` to the Cauchy transform of a
free Poisson distribution with rate `lambda>=0` and jump size `alpha` at `z`
(default `lambda=1`, `alpha=1`).
"""
function poisson(; lambda::Real = 1.0, alpha::Real = 1.0)
    if lambda < 0
        error("rate lambda has to be non-negative, but is $lambda")
    end
    function f(z::Number)
        z = complex(z)
        numroot = sqrt(z - alpha * (1 + lambda)^2 - 4 * lambda * alpha^2)
        num = z + alpha * (1 - lambda) - numroot
        return num / (2 * alpha * z)
    end
    return f
end

###############################################################################
#
#   Bernoulli distribution
#
###############################################################################

function bernoulli(p::Real)
    return bernoulli(; p)
end

"""
    bernoulli(p)
    bernoulli(; p)

Return a function mapping a complex value `z` to the Cauchy transform of a
Bernoulli distribution with probability `0<=p<=1` at `z` (default `p=0.5`).
"""
function bernoulli(; p::Real = 0.5)
    if !(0 <= p <= 1)
        error("p has to be in [0,1], but is $p")
    end
    function f(z::Number)
        z = complex(z)
        return (z - (1 - p)) / (z * (z - 1))
    end
    return f
end

###############################################################################
#
#   symmetric Bernoulli distribution
#
###############################################################################

function symmetric_bernoulli(p::Real)
    return symmetric_bernoulli(; p)
end

"""
    symmetric_bernoulli(p)
    symmetric_bernoulli(; p)

Return a function mapping a complex value `z` to the Cauchy transform of a
symmetric Bernoulli distribution with probability `0<=p<=1` at `z`
(default `p=0.5`).
"""
function symmetric_bernoulli(; p::Real = 0.5)
    if !(0 <= p <= 1)
        error("p has to be in [0,1], but is $p")
    end
    function f(z::Number)
        z = complex(z)
        return (z - (1 - 2*p)) / (z^2 - 1)
    end
    return f
end

###############################################################################
#
#   Cauchy distribution
#
###############################################################################

function cauchy(a::Real, t::Real)
    return cauchy(; a, t)
end

"""
    cauchy(a, t)
    cauchy(; a, t)

Return a function mapping a complex value `z` to the Cauchy transform of a
Cauchy distribution with location parameter `a` and scale parameter `t>0` at `z`
 (default `a=0`, `t=1`).
"""
function cauchy(; a::Real = 0.0, t::Real = 1.0)
    if t <= 0
        error("scale parameter t has to be positive, but is $t")
    end
    return z::Number -> inv(complex(z) - (a - t * im))
end

###############################################################################
#
#   arcsine distribution
#
###############################################################################

"""
    arcsine()

Return a function mapping a complex value `z` to the Cauchy transform of an
arcsine distribution at `z`.
"""
function arcsine()
    return z::Number -> inv(sqrt(complex(z)^2 - 4))
end
