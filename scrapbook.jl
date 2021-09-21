include("poly_factorization_project.jl")

function test_factor()
    Random.seed!(0)

    prime = 3
    for _ in 1:20
        @show p = PolynomialModP(rand(Polynomial), prime)
        @show factorization = factor(p)
        @show pr = expand_factorization(factorization)
        @assert iszero(p - pr)
    end
end

#test_factor()
x = x_poly()

p = PolynomialModP(x^9 + x^6 + x^3 + x, 3)

factor(p)