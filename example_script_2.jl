include("poly_factorization_project.jl")

println("This package provides a way to represent integer polynomials symbolically.")
println("Once constructed it can be used to perform arithmetic, take derivatives or factorise.")
println()

println("The basic polynomial x can be constructed as:")
@show x = x_poly()
println()

println("Arithmetic with integers is supported")
@show p1 = 3(x - 1)
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
println("Arithmetic between polynomials:")
@show p1 + p3
@show p2 - p4
@show p1 * p3


#p2_d_p1 = p2 ÷ p1
#@show p2_d_p1(11)
# mod

# division

# factor

# expand_factorization

# distinct degree