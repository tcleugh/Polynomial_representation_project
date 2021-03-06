#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

"""
A term.
"""
struct Term  #structs are immutable by default
    coeff::Int
    degree::Int
    function Term(coeff::Int, degree::Int)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

"""
Creates the zero term.
"""
zero(::Type{Term})::Term = Term(0,0)

"""
Creates the unit term.
"""
one(::Type{Term})::Term = Term(1,0)

###########
# Display #
###########

"""
Show a term.
"""
# Changes made here to outputs of a⋅x^1 → a⋅x, b⋅x^0 → b, 1⋅x^n → x^n
function show(io::IO, t::Term) 
    if t.degree == 0 
        output = "$(t.coeff)"
    elseif t.degree == 1
        output = t.coeff == 1 ? "x" : "$(t.coeff)⋅x"
    else
        output = t.coeff == 1 ? "x^$(t.degree)" : "$(t.coeff)⋅x^$(t.degree)"
    end
    print(io, output)
end

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::Term)::Bool = iszero(t.coeff)

"""
Compare two terms.
"""
isless(t1::Term,t2::Term)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  

"""
Evaluate a term at a point x.
"""
evaluate(t::Term, x::T) where T <: Number = t.coeff * x^t.degree

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1::Term,t2::Term)::Term
    @assert t1.degree == t2.degree
    Term(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::Term,) = Term(-t.coeff,t.degree)  

"""
Subtract two terms with the same degree.
"""
-(t1::Term, t2::Term)::Term = t1 + (-t2) 

"""
Multiply two terms.
"""
*(t1::Term, t2::Term)::Term = Term(t1.coeff * t2.coeff, t1.degree + t2.degree)


"""
Compute the mod of a term with an integer.
"""
mod(t::Term, p::Int)::Term = Term(mod(t.coeff,p), t.degree)

"""
Compute the symmetric mod of a term with an integer.
"""
smod(t::Term, p::Int)::Term = Term(smod(t.coeff, p), t.degree)

"""
Compute the derivative of a term.
"""
derivative(t::Term)::Term = Term(t.coeff*t.degree, max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::Term,t2::Term) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Int)::Term = Term(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::Term, n::Int) = t ÷ Term(n,0)

"""
Returns the vector sorted with all terms of same degree combined
"""
function merge(terms::Vector{Term})::Vector{Term}
    sort!(terms)
    fixed_terms = Term[]
    i = 1
    while i ≤ length(terms)
        t = terms[i]
        coeff = t.coeff
        while i < length(terms) && terms[i + 1].degree == t.degree
            coeff += terms[i + 1].coeff
            i += 1
        end
        push!(fixed_terms, Term(coeff, t.degree))
        i += 1
    end
    return fixed_terms
end
