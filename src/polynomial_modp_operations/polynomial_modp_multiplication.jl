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

    # Calculates all term multiplications
    terms = Vector{Term}(undef, length(p1) * length(p2))
    
    i = 1
    for t1 in p1
        for t2 in p2
            terms[i] = t1*t2
            i += 1
        end
    end

    # Polynomial constructor handles terms of same degree
    return PolynomialModP(terms, p1.prime)
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

    for i in 1:length(digs)
        if @inbounds digs[i] == 1 
           out *= square
        end
        (i == length(digs)) && break
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
