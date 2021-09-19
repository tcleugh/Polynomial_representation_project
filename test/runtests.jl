#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

include("../poly_factorization_project.jl")

####
# Execute unit tests for integers
###
include("integers_test.jl")
test_euclid_ints()
test_ext_euclid_ints()

####
# Execute unit tests for polynomials
####
include("polynomials_test.jl")
prod_test_poly()
prod_derivative_test_poly()
#ext_euclid_test_poly()
#division_test_poly()

####
# Execute unit tests for polynomial factorization
####
#include("factorization_test.jl")
#factor_test_poly()

####
# Execute unit tests for polynomial mod p
####
include("polynomial_modp_test.jl")
test_coeff_modp_after_construct()
test_evaluate_modp()
test_derivative_modp()
prod_test_poly_modp()
division_test_poly_modp()
ext_euclid_test_poly_modp()
power_test_poly_modp()