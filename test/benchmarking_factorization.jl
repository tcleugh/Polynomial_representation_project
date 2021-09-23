include("../poly_factorization_project.jl")


"""
Test product of polynomials.
"""
function test_factor(poly_list::Vector{PolynomialModP})

    for p in poly_list
        a = factor(p)
    end
end

function time_factor(num_primes::Int)
    Random.seed!(0)
    N = 5

    times = Vector{Float64}(undef,0)
    prime_range = Vector{Int}(undef,0)
    
    prime = 3

    for _ in 1:num_primes
        prime = nextPrime(prime)
        push!(prime_range, prime)

        poly_list = Vector{PolynomialModP}(undef, N)
        for i in 1:N
            poly_list[i] = rand(PolynomialModP, prime)
        end    

        start_time = time_ns()
        test_factor(poly_list)
        end_time = time_ns()

        push!(times, Int(end_time - start_time)/10^6) #time in milliseconds
    end
    return times, prime_range
end


function plot_timings()
    
    times, prime_range = time_factor(20)

    display(plot(prime_range , times, 
        title = "Comparing Factorization mod p", 
        xlabel = "Prime", 
        ylabel = "Factorization time (ms)", 
        shape = :circle,
        legend = false))


end

plot_timings();