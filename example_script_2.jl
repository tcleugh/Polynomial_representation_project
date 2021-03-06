include("poly_factorization_project.jl")

println("This package provides a way to represent integer polynomials symbolically.")
println("Once constructed it can be used to perform arithmetic, take derivatives or factorise.")
println()

println("The basic polynomial x can be constructed as:")
@show x = x_poly()
println()

println("Arithmetic with integers is supported")
@show p1 = x^4 + 5x^2 - 3x - 7
println()


println("A random polynomial can also be created (Note it is not uniformlly distributed)")
Random.seed!(0)
@show p2 = rand(Polynomial)
println()

println("Each polynomial is made up of a series of terms. These terms can be created seperately:")
@show t1 = one(Term)
@show t2 = Term(1, 3)
@show t3 = Term(5, 1)
println()

println("A vector of terms can be used to create a polynomial.")
@show p3 = Polynomial([t1, t2, t3])
println("New terms can also be added to a polynomial")
@show p3 += Term(2, 4)
println()

println("Some useful methods on polynomials:")
@show leading(p3)
@show coeffs(p3)
@show degree(p3)
@show evaluate(p3, 2)
@show p4 = derivative(p3)
println()

println("Polynomials support equality checking:")
@show p1 == p1
@show p1 == p2
println()

# arithmetic
println("Arithmetic with polynomials:")
@show p1 + p3
@show p2 - p4
@show p1 * p3
@show p1^3
println()

println("Polynomials can also be created in a field Z mod prime")
@show pm1 = PolynomialModP(23x^2 +7x + 8, 7)
@show pm2 = rand(PolynomialModP, 7)
println()

println("Modular arithmetic can be performed on polynomials in the same field")
@show pm1 + pm2
@show pm1 - pm2
@show pm1 * pm2
println()

println("Polynomials mod p also support divison and gcd")
@show gcd(pm1, pm2)
@show (pm2 ÷ pm1)
println()

println("Polynomials mod p can also be factored mod p")
@show f = factor(pm2)
@show expand_factorization(f);
