#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials.
"""
function *(p1::Polynomial, p2::Polynomial)::Polynomial
    p_out = Polynomial()
    for t in p1
        p_out = p_out + (t * p2)
    end
    return p_out
end

"""
Power of a polynomial
"""
function ^(p::Polynomial, n::Integer)
    n < 0 && error("No negative power")
    out = one(p)
    n == 0 && return out
    
    digs = digits(n, base=2)
    len = length(digs)
    square = p

    for i in 1:(len)
        if digs[i] == 1 
           out *= square
        end
        (i == 0) && break
        square *= square
    end
    return out
end


"""
Multiply two polynomials using chinese remainder therom
"""
function crt_mult(a::Polynomial, b::Polynomial)::Polynomial
    (iszero(a) || iszero(b)) && return zero(Polynomial)
    #(degree(a) == 0 || degree(b) == 0) && return a * b

    h1, h2 = maximum(abs.(coeffs(a))), maximum(abs.(coeffs(b)))

    B = 2 * h1 * h2 * min(degree(a) + 1, degree(b) + 1)
    M = 3
    c = PolynomialModP(a, 3) * PolynomialModP(b, 3)
    prime = M

    while M < B
        prime = nextPrime(prime)
        c_ = PolynomialModP(a, prime) * PolynomialModP(b, prime)
        c = crt(c, c_, M, prime)
        M = M * prime
    end
    
    return smod(c, M)
end

"""
Chinese remainder therom on two polynomials
"""
function crt(p1::Polynomial, p2::Polynomial, n::Integer, m::Integer)::Polynomial
    d1, d2 = degree(p1), degree(p2)

    out = zero(Polynomial)
    
    k = max(d1, d2)
    while k ≥ 0
        ak = (k > d1) ? 0 : coeff(p1, k)
        bk = (k > d2) ? 0 : coeff(p2, k)
        ck = crt(ak, bk, n, m)
        push!(out, Term(ck, k))
        k -= 1
    end
    return out
end

crt(p1::PolynomialModP, p2::PolynomialModP, n::Integer, m::Integer)::Polynomial = crt(p1.poly, p2.poly, n, m)
crt(p1::Polynomial, p2::PolynomialModP, n::Integer, m::Integer)::Polynomial = crt(p1, p2.poly, n, m)
crt(p1::PolynomialModP, p2::Polynomial, n::Integer, m::Integer)::Polynomial = crt(p1.poly, p2, n, m)

"""
Multiplication of polynomial and term.
"""
*(t::Term, p::Polynomial)::Polynomial = iszero(t) ? Polynomial() : Polynomial(sort(map((term) -> t * term, p.terms)))
*(p::Polynomial, t::Term)::Polynomial = t * p

"""
Multiplication of polynomial and an integer.
"""
*(n::Int, p::Polynomial)::Polynomial = p * Term(n, 0)
*(p::Polynomial, n::Int)::Polynomial = n * p