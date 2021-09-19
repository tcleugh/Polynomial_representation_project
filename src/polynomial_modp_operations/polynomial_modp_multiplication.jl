#############################################################################
#############################################################################
#
# This file implements polynomial mod p Multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials mod p.
"""
function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    p_out = zero(PolynomialModP, p1.prime)

    for t in p1
        p_out += PolynomialModP(t * p2.poly, p1.prime)
    end
    return p_out
end

"""
Power of a polynomial mod p.
"""
function ^(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
    out = one(PolynomialModP, p.prime)
    n == 0 && return out
    
    digs = digits(n, base=2)
    len = length(digs)
    square = p

    for i in 1:(len)
        if digs[i] == 1 
           out *= square
        end
        (i == 0) && break
        square *= square
    end
    return out
end

"""
(Outdated: only use for benchmarking new alg) Power of a polynomial mod p.
"""
function old_pow(p::PolynomialModP, n::Int)
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
