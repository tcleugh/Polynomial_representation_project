include("../poly_factorization_project.jl")


"""
Test product of polynomials.
"""
function test_regular_mult(degree::Int, max_coeff::Int; N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(Polynomial, degree = degree, max_coeff = max_coeff)
        p2 = rand(Polynomial, degree = degree, max_coeff = max_coeff)
        prod = p1*p2
    end
end

"""
Test product of polynomials using crt multiplication.
"""
function test_crt_mult(degree::Int, max_coeff::Int; N::Int = 10^3, N_prods::Int = 10, seed::Int = 0)
    
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(Polynomial, degree = degree, max_coeff = max_coeff)
        p2 = rand(Polynomial, degree = degree, max_coeff = max_coeff)
        prod = crt_mult(p1, p2)
    end
end

function time_regular_vary_degree(max_degree::Int)
    
    times = Vector{Float64}(undef,0)
    degree_range = 5: 50 :max_degree

    for degree in degree_range
        start_time = time_ns()
        test_regular_mult(degree, 100)
        end_time = time_ns()
        push!(times, Int(end_time - start_time)/10^6) #time in milliseconds
    end
    return times, degree_range
end

function time_crt_vary_degree(max_degree::Int)
    
    times = Vector{Float64}(undef,0)
    degree_range = 5: 50 : max_degree

    for degree in degree_range
        start_time = time_ns()
        test_crt_mult(degree, 100)
        end_time = time_ns()
        push!(times, Int(end_time - start_time)/10^6) #time in milliseconds
    end
    return times, degree_range
end

function time_regular_vary_coeff(max_power::Int)
    
    times = Vector{Float64}(undef,0)
    coeff_range = [5^x for x in 1:max_power]

    for coeff in coeff_range
        start_time = time_ns()
        test_regular_mult(10, coeff)
        end_time = time_ns()
        push!(times, Int(end_time - start_time)/10^6) #time in milliseconds
    end
    return times, coeff_range
end

function time_crt_vary_coeff(max_power::Int)
    
    times = Vector{Float64}(undef,0)
    coeff_range = [5^x for x in 1:max_power]

    for coeff in coeff_range
        start_time = time_ns()
        test_crt_mult(10, coeff)
        end_time = time_ns()
        push!(times, Int(end_time - start_time)/10^6) #time in milliseconds
    end
    return times, coeff_range
end


function plot_timings()
    
    times1, degree_range1 = time_regular_vary_degree(1000)
    times2, degree_range2 = time_crt_vary_degree(1000)

    display(plot([degree_range1 degree_range2], [times1 times2], 
        title = "Comparing methods with varying degree",
        label = ["Regular" "CRT"], 
        xlabel = "Degree of polynomials", 
        ylabel = "Multiplication time (ms)", 
        shape = :circle,
        legend = :topleft))


    times3, coeff_range1 = time_regular_vary_coeff(10)
    times4, coeff_range2 = time_crt_vary_coeff(10)

    display(plot([coeff_range1 coeff_range2], [times3 times4], 
        title = "Comparing with varying maximum coefficient",
        label = ["Regular" "CRT"], 
        xlabel = "Maximum coefficients", 
        ylabel = "Multiplication time (ms)", 
        shape = :circle,
        legend = :topleft,
        xscale = :log10))    

end

plot_timings();