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

println()
println()

@show g1,s1,t1 = extended_euclid_alg(p1, p2, 11)
@show mod(s1*p1 + t1*p2 - g1, 11)

println()
p1m = PolynomialModP(p1, 11)
p2m = PolynomialModP(p2, 11)

@show g,s,t = extended_euclid_alg(p1m, p2m)

@show mod(s*p1m + t*p2m - g, 11)
