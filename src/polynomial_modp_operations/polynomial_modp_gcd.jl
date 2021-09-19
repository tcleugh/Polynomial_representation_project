#############################################################################
#############################################################################
#
# This file implements polynomial mod p GCD operations
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials mod p
"""
function extended_euclid_alg(p1::PolynomialModP, p2::PolynomialModP)::Tuple{PolynomialModP, PolynomialModP, PolynomialModP}
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
gcd(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP = extended_euclid_alg(p1,p2) |> first

"""
Take the mod of a polynomial mod p with an integer.
"""
mod(p::PolynomialModP, n::Int)::PolynomialModP = PolynomialModP(mod(p.poly, n), p.prime)
