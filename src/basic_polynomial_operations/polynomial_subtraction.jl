#############################################################################
#############################################################################
#
# This file implements polynomial subtraction 
#                                                                               
#############################################################################
#############################################################################

"""
Subtraction of two polynomials.
"""
-(p1::Polynomial, p2::Polynomial)::Polynomial = p1 + (-p2)

"""
Subtraction of polynomial and a term
"""
-(p::Polynomial, t::Term)::Polynomial = p - Polynomial(t)
-(t::Term, p::Polynomial)::Polynomial = Polynomial(t) - p

"""
Subtraction of polynomial and an integer
"""
-(p::Polynomial, n::Int)::Polynomial = p - Term(n,0)
-(n::Int, p::Polynomial)::Polynomial = Term(n,0) - p 