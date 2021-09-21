#############################################################################
#############################################################################
#
# This file implements polynomial mod p factorization 
#                                                                               
#############################################################################
#############################################################################

"""
Factors a polynomial over the field Z_p.

Returns a vector of tuples of (irreducible polynomials (mod p), multiplicity) such that their product of the list (mod p) is f. Irreducibles are fixed points on the function factor.
"""
function factor(p::PolynomialModP)::Vector{Tuple{PolynomialModP,Integer}}
    #Cantor Zassenhaus factorization
    degree(p) ≤ 1 && return [(p, 1)]

    # make p primitive
    pp = prim_part(p)      
    #@show "after prim:", pp

    # make p square-free
    squares_poly = gcd(p, derivative(pp)) 
    pp = pp ÷ squares_poly 
    #@show "after square free:", pp

    # make p monic
    pp = pp ÷ leading(pp).coeff     
    #@show "after monic:", pp

    dds = dd_factor(pp)

    ret_val = Tuple{PolynomialModP,Int}[]

    for (k, dd) in enumerate(dds)
        sp = dd_split(dd, k)
        sp = map((p) -> p ÷ leading(p).coeff ,sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, multiplicity(p, mp)) )
        end
    end

    #Append the leading coefficient as well
    push!(ret_val, (PolynomialModP(Term(leading(p).coeff, 0), p.prime), 1))

    return ret_val
end

"""
Expand a factorization.
"""
function expand_factorization(factorization::Vector{Tuple{PolynomialModP, Integer}})::PolynomialModP
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

"""
Compute the number of times g divides f
"""
function multiplicity(f::PolynomialModP, g::PolynomialModP)::Integer
    degree(gcd(f, g)) == 0 && return 0
    return 1 + multiplicity(f ÷ g, g)
end


"""
Distinct degree factorization.

Given a square free polynomial `p` returns a list, `g` such that `g[k]` is a product of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the product of the list (mod `prime`) is equal to `p` (mod `prime`).
"""
function dd_factor(p::PolynomialModP)::Array{PolynomialModP}
    x = PolynomialModP(x_poly(), p.prime)
    w = x
    g = Array{PolynomialModP}(undef, degree(p)) #Array of polynomials indexed by degree

    #Looping over degrees
    for k in 1:degree(p)
        w = rem(w^p.prime, p) # This line causes issues for mod 3
        g[k] = gcd(w - x, p) 
        p = p ÷ g[k]
    end

    #edge case for final factor
    p != one(PolynomialModP, p.prime) && push!(g, p)
    return g

end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of that list (mod prime) is the polynomial `f`.
"""
function dd_split(p::PolynomialModP, d::Integer)::Vector{PolynomialModP}
    degree(p) == d && return [p]
    degree(p) == 0 && return []
    w = PolynomialModP(rand(Polynomial, degree = d, monic = true), p.prime)
    n_power = (p.prime^d - 1) ÷ 2
    g = gcd(w^n_power - one(PolynomialModP, p.prime), p)
    ḡ = p ÷ g 
    return vcat(dd_split(g, d), dd_split(ḡ, d))
end

"""
Pretty printing of the factored polynomials
"""
function show(io::IO, factoredp::Vector{Tuple{PolynomialModP,Integer}})
    length(factoredp) == 0 && print(io, 0)
    
    print(io, first(last(factoredp)).poly)

    for i in 1:(length(factoredp) - 1)
        (p, n) = factoredp[i]
        
        print(io, "($(p.poly))")
        n > 1 && print(io, "^$n")
    end

    print(io, " mod $(first(factoredp[1]).prime)")
end