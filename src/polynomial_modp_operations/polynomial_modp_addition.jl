#############################################################################
#############################################################################
#
# This file implements polynomial mod p addition 
#                                                                               
#############################################################################
#############################################################################


"""
Addition of two polynomials mod p.
"""
function +(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    return PolynomialModP(p1.poly + p2.poly, p1.prime)
end 

"""
Add a polynomial and a term mod p.
"""
+(p::PolynomialModP, t::Term) = PolynomialModP(p.poly + t, p.prime)
+(t::Term, p::PolynomialModP) = p + t

"""
Add a polynomial and an integer.
"""
+(p::PolynomialModP, n::Int) = PolynomialModP(p.poly + n, p.prime)
+(n::Int, p::PolynomialModP) = n + p