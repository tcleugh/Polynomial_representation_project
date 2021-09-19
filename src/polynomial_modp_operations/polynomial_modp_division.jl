#############################################################################
#############################################################################
#
# This file implements polynomial mod p division 
#                                                                               
#############################################################################
#############################################################################

"""
Integer division of two polynomials mod p.
"""
function divide(num::PolynomialModP, den::PolynomialModP)
    @assert num.prime == den.prime
    rem, prime = num, num.prime

    degree(rem) < degree(den) && return (zero(PolynomialModP, prime), rem)
    iszero(den) && throw(DivideError())
    quo = PolynomialModP(prime)
    prev_degree = degree(rem)
    while degree(rem) ≥ degree(den) 
        h = PolynomialModP((leading(rem) ÷ leading(den))(prime), prime)
        rem = rem - h * den
        quo = quo + h  
        prev_degree == degree(rem) && break
        prev_degree = degree(rem)
    end
    @assert iszero(num  - (quo * den + rem))
    return quo, rem

end

"""
The quotient from polynomial division mod p.
"""
÷(num::PolynomialModP, den::PolynomialModP)  = first(divide(num, den))

"""
The remainder from polynomial division mod p.
"""
rem(num::PolynomialModP, den::PolynomialModP)  = last(divide(num,den))

"""
Integer division of a polynomial by an integer mod p.
"""
÷(p::PolynomialModP, n::Int)::PolynomialModP = p ÷ PolynomialModP(Term(n, 0), p.prime)  # PolynomialModP(p.poly ÷ n, p.prime)