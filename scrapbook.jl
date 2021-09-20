include("poly_factorization_project.jl")

x = x_poly()

p = PolynomialModP(x^2 + 2x + 1, 11)

@show factor(p)