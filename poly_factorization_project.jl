#############################################################################
#############################################################################
#
# This is the main project file for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

using DataStructures, Distributions, StatsBase, Random, BenchmarkTools

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one

include("src/general_alg.jl")
include("src/term.jl")
include("src/polynomial.jl")
include("src/polynomial_mod_p.jl")

    include("src/basic_polynomial_operations/polynomial_addition.jl")
    include("src/basic_polynomial_operations/polynomial_subtraction.jl")
    include("src/basic_polynomial_operations/polynomial_multiplication.jl")

    include("src/polynomial_modp_operations/polynomial_modp_addition.jl")
    include("src/polynomial_modp_operations/polynomial_modp_subtraction.jl")
    include("src/polynomial_modp_operations/polynomial_modp_multiplication.jl")
    include("src/polynomial_modp_operations/polynomial_modp_factorization.jl")
    include("src/polynomial_modp_operations/polynomial_modp_division.jl")
    include("src/polynomial_modp_operations/polynomial_modp_gcd.jl")