#############################################################################
#############################################################################
#
# This file implements polynomial mod p Subtraction 
#                                                                               
#############################################################################
#############################################################################

"""
Subtraction of two polynomials mod p.
"""
-(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP = p1 + (-p2)
 
"""
Subtraction of polynomial and a term mod p
"""
-(p::PolynomialModP, t::Term)::PolynomialModP = PolnomialModP(p.poly - t, p.prime)
-(t::Term, p::PolynomialModP)::PolynomialModP = PolnomialModP(t - p.poly, p.prime)

"""
Subtraction of polynomial and an integer mod p
"""
-(p::PolynomialModP, n::Int) = PolnomialModP(p.poly - n, p.prime)
-(n::Int, p::PolynomialModP) = PolnomialModP(n - p.poly, p.prime)
