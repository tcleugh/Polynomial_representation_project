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
function PolynomialModP(term::Term, prime::Integer)
    return PolynomialModP(Polynomial(term), prime)
end

"""
Construct the zero polynomial in the given field.
"""
function PolynomialModP(prime::Integer)
    return PolynomialModP(Polynomial(), prime)
end

"""
Construct the zero polynomial in the given field.
"""
function zero(::Type{PolynomialModP}, prime::Integer)
    return PolynomialModP(prime)
end

"""
Construct the one polynomial in the given field.
"""
function one(::Type{PolynomialModP}, prime::Integer)
    return PolynomialModP(one(Polynomial), prime)
end

"""
Generates a random polynomial mod p.
"""
function rand(::Type{PolynomialModP}, prime::Integer)
    return PolynomialModP(rand(Polynomial), prime)
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
derivative(p::PolynomialModP) = PolynomialModP(derivative(p.poly), p.prime)

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialModP) = PolynomialModP(p.poly ÷ content(p), p.prime)

"""
A square free polynomial.
"""
square_free(p::PolynomialModP)::PolynomialModP = PolnomialModP(square_free(p.poly, p.prime), p.prime)


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

#####################################
# Operations with polynomials mod p #
#####################################

### Addition ###

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

### Subtraction ###

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

### Multiplication ###

"""
Multiply two polynomials mod p.
"""
function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    return PolynomialModP(p1.poly * p2.poly, p1.prime)
end

"""
Multiply two polynomials mod p.
"""
function old_mult(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    return PolynomialModP(p1.poly * p2.poly, p1.prime)
end

"""
Power of a polynomial mod p.
"""
function ^(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
    out = one(PolynomialModP, p.prime)
    for _ in 1:n
        out*= p
    end
    return out
end

"""
Multiplication of polynomial and term.
"""
*(t::Term, p::PolynomialModP)::PolynomialModP = PolynomialModP(t * p.poly, p.prime)
*(p::PolynomialModP, t::Term)::PolynomialModP = t * p

"""
Multiplication of polynomial and an integer.
"""
*(n::Int, p::PolynomialModP)::PolynomialModP = PolynomialModP(p.poly * n, p.prime)
*(p::PolynomialModP, n::Int)::PolynomialModP = n * p

### Division ###

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
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialModP, n::Int)::PolynomialModP = PolynomialModP(p.poly ÷ n, p.prime)

"""
The extended euclid algorithm for polynomials mod p
"""
function extended_euclid_alg1(p1::PolynomialModP, p2::PolynomialModP)
    @assert p1.prime == p2.prime
    a, b, prime = p1.poly, p2.poly, p1.prime

    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(Polynomial), zero(Polynomial)
    old_t, t = zero(Polynomial), one(Polynomial)

    while !iszero(mod(r, prime))
        q = divide(old_r, r)(prime) |> first
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t
    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

"""
The extended euclid algorithm for polynomials mod p
"""
function extended_euclid_alg(p1::PolynomialModP, p2::PolynomialModP)
    @assert p1.prime == p2.prime
    prime = p1.prime

    old_r, r = p1, p2
    old_s, s = one(PolynomialModP, prime), zero(PolynomialModP, prime)
    old_t, t = zero(PolynomialModP, prime), one(PolynomialModP, prime)

    while !iszero(r)
        q = first(divide(old_r, r)) 
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    end
    g, s, t = old_r, old_s, old_t
    @assert iszero(s*p1 + t*p2 - g)
    return g, s, t  
end

"""
The GCD of two polynomials modulo prime.
"""
gcd(p1::PolynomialModP, p2::PolynomialModP) = extended_euclid_alg(p1,p2) |> first

"""
Take the mod of a polynomial mod p with an integer.
"""
mod(p::PolynomialModP, n::Int)::PolynomialModP = PolynomialModP(mod(p.poly, n), p.prime)
