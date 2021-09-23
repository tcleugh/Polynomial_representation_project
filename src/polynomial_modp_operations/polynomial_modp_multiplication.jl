#############################################################################
#############################################################################
#
# This file implements polynomial mod p Multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials mod p
"""
function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    (iszero(p1) || iszero(p2)) && return zero(PolynomialModP, p1.prime)

    prime = p1.prime
    max_degree = degree(p1) + degree(p2)

    # Creates a sorted vector of the coefficients of p1 * p2
    coeffs = fill(0, (max_degree + 1, 1))
    for t1 in p1
        for t2 in p2
            @inbounds coeffs[t1.degree + t2.degree + 1] += t1.coeff * t2.coeff
        end
    end

    # Converts the coefficient vector into a term vector
    fixed_terms = Term[]
    for (degree, coeff) in enumerate(coeffs)
        coeff != 0 && mod(coeff, prime) != 0 && push!(fixed_terms, Term(mod(coeff, prime), degree - 1))
    end
    # uses the "safe" constructor to avoid excess sorting/merging/filtering
    return PolynomialModP(fixed_terms, p1.prime, true)
end

"""
Power of a polynomial mod p.
"""
function ^(p::PolynomialModP, n::Int)::PolynomialModP
    n < 0 && error("No negative power")
    out = one(PolynomialModP, p.prime)
    n == 0 && return out
    
    digs = digits(n, base=2)
    square = p
    len = length(digs)

    for i in 1:len
        if @inbounds digs[i] == 1 
           out *= square
        end
        (i == len) && break
        square *= square
    end
    return out
end

"""
(Outdated: only use for benchmarking new alg) Power of a polynomial mod p.
"""
function old_pow(p::PolynomialModP, n::Int)::PolynomialModP
    n < 0 && error("No negative power")
    out = one(PolynomialModP, p.prime)
    for _ in 1:n
        out*= p
    end
    return out
end

"""
Multiplication of polynomial and term mod p.
"""
*(t::Term, p::PolynomialModP)::PolynomialModP = PolynomialModP(t * p.poly, p.prime)
*(p::PolynomialModP, t::Term)::PolynomialModP = t * p

"""
Multiplication of polynomial and an integer mod p.
"""
*(n::Int, p::PolynomialModP)::PolynomialModP = PolynomialModP(p.poly * n, p.prime)
*(p::PolynomialModP, n::Int)::PolynomialModP = n * p
