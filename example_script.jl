include("poly_factorization_project.jl")

x = x_poly()
p1 = 2x^3 + 4x^2 - 3x
p2 = 2x^4 - 4x^2 - 3x + 3

@show p1*p2
@show p1^3

@show derivative(p1*p2)
@show derivative(p1)*p2 + p1*derivative(p2);

prime = 17
p = mod((7x^3 + 2x^2 + 8x + 1)*(x^2+x+1),prime)
println("Will factor this polynomial (mod $prime): ", p)
factorization = factor(p,prime)
println("Here is the factorization: ", factorization)

pr = mod(expand_factorization(factorization),prime)
println("Reconstructing: ", pr)

@show p4 = PolynomialModP(p1, 3)