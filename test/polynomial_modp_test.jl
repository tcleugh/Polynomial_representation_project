#############################################################################
#############################################################################
#
# This file contains units tests for polynomial mod p operations
#                                                                               
#############################################################################
#############################################################################

"""
Test creation of polynomial mod p
"""
function test_coeff_modp_after_construct(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        for i in [3, 5, 7, 11, 17]
            p = rand(PolynomialModP, i)
            for coeff in coeffs(p)
                @assert coeff != 0
                @assert mod(coeff, p.prime) == coeff
            end
        end
    end

    println("test_coeff_modp_after_construct - PASSED")
end

"""
Test evaluation of polynomial mod p
"""
function test_evaluate_modp(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        for i in [3, 5, 7, 11, 17]
            p = rand(PolynomialModP, i)
            for x in -5:5
                length(p) != 0 && @assert evaluate(p, x) == mod(evaluate(p.poly, x), p.prime)
            end
        end
    end

    println("test_evaluate_modp - PASSED")
end

"""
Test evaluation of polynomial mod p
"""
function test_derivative_modp(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        for i in [3, 5, 7, 11, 17]
            p = rand(PolynomialModP, i)
            @assert derivative(p).poly == mod(derivative(p.poly), p.prime)
        end
    end

    println("test_derivative_modp - PASSED")
end

"""
Test product of polynomials mod p.
"""
function prod_test_poly_modp(;N::Int = 10^2, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        for prime in [3, 5, 7, 11, 17]
            p1 = rand(PolynomialModP, prime)
            p2 = rand(PolynomialModP, prime)
            prod = p1 * p2
            @assert leading(prod) == mod(leading(p1) * leading(p2), prime)
        end
    end

    for _ in 1:N
        for prime in [3, 5, 7, 11, 17]
            p_base = PolynomialModP(Term(1,0), prime)
            for _ in 1:N_prods
                p = rand(PolynomialModP, prime)
                prod = p_base * p
                @assert leading(prod) == mod(leading(p_base)*leading(p), prime)
                p_base = prod
            end
        end
    end
    println("prod_test_poly_modp - PASSED")
end

"""
Test division of polynomials mod p.
"""
function division_test_poly_modp(;prime::Int = 101, N::Int = 10^4, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        p_prod = p1*p2
        q, r = PolynomialModP(prime), PolynomialModP(prime)
        try
            q, r = divide(p_prod, p2)
            if iszero(q)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero(mod(q*p2+r - p_prod, prime))
    end
    println("division_test_poly_modp - PASSED")
end

"""
Test the extended euclid algorithm for polynomials mod p.
"""
function ext_euclid_test_poly_modp(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        g, s, t = extended_euclid_alg(p1, p2)
        @assert iszero(s*p1 + t*p2 - g)
    end
    println("ext_euclid_test_poly_modp - PASSED")
end


"""
Test the power function for polynomials mod p.
"""
function power_test_poly_modp(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialModP, prime)
        
        n = rand(0:20)
        pn = one(PolynomialModP, prime)
        for _ in 1:n
            pn *= p
        end

        @assert p^n == pn
    end
    println("power_test_poly_modp - PASSED")
end

