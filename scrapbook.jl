include("poly_factorization_project.jl")

x = x_poly()

#Random.seed!(0)

@show p1 = rand(Polynomial)
@show p2 = rand(Polynomial)

@show p1*p2
