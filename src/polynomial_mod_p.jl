#############################################################################
#############################################################################
#
# This file defines the polynomial mod P type with several operations 
#                                                                               
#############################################################################
#############################################################################

##########################################
# Polynomial mod p type and construction #
##########################################

"""
A Polynomial type - designed to be for polynomials with integer coefficients.
"""
struct PolynomialModP
    poly::Polynomial   
    prime::Integer
    PolynomialModP(f, p) = new(Polynomial(map((term) -> mod(term, p), f.terms)) , p)
end

"""
Construct the zero polynomial in the given field.
"""
function PolynomialModP(prime::Integer)
    return PolynomialModP(Polynomial(), prime)
end

###########
# Display #
###########

"""
Show a polynomial mod P.
"""
function show(io::IO, p::PolynomialModP) 
    print(io, p.poly)
    print(io, " mod $(p.prime)")
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the terms of the polynomial. The iteration is in an arbitrary order.
"""
iterate(p::PolynomialModP, state=1) = iterate(p.poly, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialModP) = length(p.poly)

"""
The leading term of the polynomial.
"""
leading(p::PolynomialModP)::Term = leading(p.poly)

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::PolynomialModP)::Vector{Int} = coeffs(p.poly)

"""
The degree of the polynomial.
"""
degree(p::PolynomialModP)::Int = degree(p.poly)

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialModP)::Int = content(p.poly)

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::PolynomialModP, x::T) where T <: Number = evaluate(p.poly, x)
