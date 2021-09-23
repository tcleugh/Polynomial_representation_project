#############################################################################
#############################################################################
#
# This file defines the polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial type and construction #
####################################

"""
A Polynomial type - designed to be for polynomials with integer coefficients.
"""
struct Polynomial
    terms::Vector{Term}   
        #The terms in the vector need to satisfy:
            # Will never have terms with 0 coefficient
            # Will never have two terms with same degree
            # Will always be sorted 
    Polynomial() = new(Term[])

    #Inner constructor
    Polynomial(h::Vector{Term}) = new(merge(filter((term) -> term.coeff != 0, h)))

    # Construct a polynomial from a vector of terms that is pre sorted with all valid coefficients
    # Only use from multiplication (allows for speed increase)
    Polynomial(h::Vector{Term}, safe::Bool) = safe ? new(h) : error("Invalid use")
end

"""
Construct a polynomial with a single term.
"""
function Polynomial(t::Term)::Polynomial
    terms = Term[]
    t.coeff != 0 && push!(terms, t)
    return Polynomial(terms)
end

"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomial(p::Int)::Polynomial = Polynomial([Term(1,p), Term(-1,0)])

"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial(n::Int)::Polynomial = Polynomial([Term(1,1), Term(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_poly()::Polynomial = Polynomial(Term(1,1))

"""
Creates the zero polynomial.
"""
zero(::Type{Polynomial})::Polynomial = Polynomial()

"""
Creates the unit polynomial.
"""
one(::Type{Polynomial})::Polynomial = Polynomial(one(Term))
one(p::Polynomial) = one(typeof(p))

"""
Generates a random polynomial.
"""
function rand(::Type{Polynomial} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree, prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1, _terms, replace = false)), _degree)
        coeffs = rand(1:max_coeff, _terms+1)
        monic && (coeffs[end] = 1)
        p = Polynomial( [Term(coeffs[i], degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::Polynomial) 
    p = deepcopy(p)
    if iszero(p)
        print(io, "0")
    else
        n = length(p.terms)
        print(io, p.terms[n])
        n -= 1
        while n > 0
            print(io, p.terms[n].coeff >= 0 ? " + " : " ", p.terms[n])
            n -= 1
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the terms of the polynomial. The iteration is in an arbitrary order.
"""
iterate(p::Polynomial, state=1) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::Polynomial)::Integer = length(p.terms)

"""
The leading term of the polynomial.
"""
leading(p::Polynomial)::Term = isempty(p.terms) ? zero(Term) : last(p.terms) 

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::Polynomial)::Vector{Int} = [t.coeff for t in p]

"""
Returns the coefficient of the term in the polynomial with given degree
"""
function coeff(p::Polynomial, k::Integer)::Integer
    for term in p.terms
        if term.degree == k
            return term.coeff
        end
    end
    return 0
end

"""
The degree of the polynomial.
"""
degree(p::Polynomial)::Integer = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::Polynomial)::Integer = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(p::Polynomial, x::T) where T <: Number = sum(evaluate(term, x) for term in p)

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
function push!(p::Polynomial, t::Term) 
    iszero(t) && return #don't push a zero
    insorted(t, p.terms) && @assert ArgumentError("Polynomial can't have two terms of the same degree")
    # inserts term into appropriate position
    insert!(p.terms, first(searchsorted(p.terms, t)), t)
    return p
end

"""
Pop the leading term out of the polynomial.
"""
pop!(p::Polynomial)::Term = pop!(p.terms)

"""
Check if the polynomial is zero.
"""
iszero(p::Polynomial)::Bool = isempty(p.terms)

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::Polynomial) = Polynomial(map((term) -> -term, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::Polynomial)::Polynomial 
    der_p = Polynomial()
    for term in p
        push!(der_p, derivative(term)) #if coeff reduced to zero, push! will handle it
    end
    return der_p
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::Polynomial) = p ÷ content(p)


#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
==(p1::Polynomial, p2::Polynomial)::Bool = p1.terms == p2.terms

"""
Check if a polynomial is equal to 0.
"""
# changed to resolve equality with integers
==(p::Polynomial, n::T) where T <: Real = (n == 0) ? iszero(p) == iszero(n) : (n == n ÷ 1) && p == Polynomial(Term(n))

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::Polynomial, n::Int) = (prime) -> Polynomial(sort(map((term) -> ((term ÷ n)(prime)), p.terms)))

"""
Take the mod of a polynomial with an integer.
"""
function mod(p::Polynomial, n::Int)::Polynomial
    p_out = Polynomial()
    for term in p
        push!(p_out, mod(term, n)) #if coeff reduced to zero, push! will handle it
    end
    return p_out
end

"""
Take the symmetric mod of a polynomial with an integer.
"""
function smod(p::Polynomial, n::Int)::Polynomial
    p_out = Polynomial()
    for term in p
        push!(p_out, smod(term, n)) #if coeff reduced to zero, push! will handle it
    end
    return p_out
end