include("poly_factorization_project.jl")

x = x_poly()

p1 = 3x - 4
p2 = 6x + 5

@show p1 * p2
@show new_mult(p1, p2)