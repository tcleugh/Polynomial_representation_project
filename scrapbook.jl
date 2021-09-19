include("poly_factorization_project.jl")

x = x_poly()

p = rand(PolynomialModP, 101)
println(p)
println(p^123)