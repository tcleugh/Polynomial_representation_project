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
    PolynomialModP(p::Polynomial, prime::Integer) = new(Polynomial(map((term) -> mod(term, prime), p.terms)), prime)
end

"""
Construct a polynomial from a single term in the given field.
"""
PolynomialModP(term::Term, prime::Integer) = PolynomialModP(Polynomial(term), prime)

"""
Construct the zero polynomial in the given field.
"""
PolynomialModP(prime::Integer) = PolynomialModP(Polynomial(), prime)

"""
Construct a polynomial from a vector of terms in the given field.
"""
PolynomialModP(h::Vector{Term}, prime::Integer) = PolynomialModP(Polynomial(h), prime)

"""
Construct the zero polynomial in the given field.
"""
zero(::Type{PolynomialModP}, prime::Integer) = PolynomialModP(prime)

"""
Construct the one polynomial in the given field.
"""
one(::Type{PolynomialModP}, prime::Integer) = PolynomialModP(one(Polynomial), prime)

"""
Generates a random polynomial mod p.
"""
rand(::Type{PolynomialModP}, prime::Integer) = PolynomialModP(rand(Polynomial), prime)


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

####################################
# Queries about a polynomial mod p #
####################################

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
Returns the coefficient mod p of the term in the polynomial with given degree
"""
coeff(p::PolynomialModP, k::Integer)::Integer = coeff(p.poly, k)

"""
The degree of the polynomial.
"""
degree(p::PolynomialModP)::Integer = degree(p.poly)

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialModP)::Integer = content(p.poly)

"""
Evaluate the polynomial at a point `x`.
"""
function evaluate(p::PolynomialModP, x::T) where T <: Number 
    total = 0
    for term in p.poly.terms
        total = mod(total + evaluate(term, x), p.prime)
    end
    return total
end

################################
# Pushing and popping of terms #
################################

"""
Push a new term mod p into the polynomial.
"""
function push!(p::PolynomialModP, t::Term) 
    iszero(t) && return #don't push a zero
    t = mod(t, p.prime)
    iszero(t) && return #don't push a zero
    push!(p.poly, t)
    return p
end

"""
Pop the leading term out of the polynomial.
"""
pop!(p::PolynomialModP)::Term = pop!(p.poly)

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialModP)::Bool = iszero(p.poly)


#############################################################################
# Transformation of the polynomial mod p to create another polynomial mod p #
#############################################################################

"""
The negative of a polynomial mod p.
"""
-(p::PolynomialModP) = PolynomialModP(-p.poly, p.prime)

"""
Create a new polynomial mod p which is the derivative of the polynomial mod p.
"""
derivative(p::PolynomialModP)::PolynomialModP = PolynomialModP(derivative(p.poly), p.prime)

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialModP)::PolynomialModP = p ÷ content(p)

"""
A square free polynomial.
"""
square_free(p::PolynomialModP)::PolynomialModP = p ÷ gcd(p, derivative(p))

#######################################
# Queries about two polynomials mod p #
#######################################

"""
Check if two polynomials are the same mod p
"""
==(p1::PolynomialModP, p2::PolynomialModP)::Bool = p1.prime == p2.prime && p1.poly == p2.poly

"""
Check if a polynomial is equal to an integer mod p.
"""
==(p::PolynomialModP, n::T) where T <: Real = (n == 0) ? iszero(p) : (n == n ÷ 1) && p == PolynomialModP(Term(n), p.prime)
