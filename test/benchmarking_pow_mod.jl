include("../poly_factorization_project.jl")


function test_new_pow(power::Int, poly_list::Vector{PolynomialModP})
    for p in poly_list
        a = p^power
    end
end

function test_old_pow(power::Int, poly_list::Vector{PolynomialModP})
    for p in poly_list
        a = old_pow(p, power)
    end
end

function time_new_vary_power(max_power::Int, poly_list::Vector{PolynomialModP})
    
    times = Vector{Float64}(undef,0)
    power_range = 0: 20 :max_power

    for power in power_range
        start_time = time_ns()
        test_new_pow(power, poly_list)
        end_time = time_ns()
        push!(times, Int(end_time - start_time)/10^6) #time in milliseconds
    end

    return times, power_range
end

function time_old_vary_power(max_power::Int, poly_list::Vector{PolynomialModP})
    
    times = Vector{Float64}(undef,0)
    power_range = 0: 20 :max_power

    for power in power_range
        start_time = time_ns()
        test_old_pow(power, poly_list)
        end_time = time_ns()
        push!(times, Int(end_time - start_time)/10^6) #time in milliseconds
    end

    return times, power_range
end


function plot_pow_timings()
    Random.seed!(0)
    N = 10^2
    max_pow = 100
    poly_list = Vector{PolynomialModP}(undef, N)
    for i in 1:N
        poly_list[i] = rand(PolynomialModP, 101)
    end

    test_new_pow(1, poly_list)
    test_old_pow(1, poly_list)

    times1, power_range1 = time_new_vary_power(max_pow, poly_list)
    times2, power_range2 = time_old_vary_power(max_pow, poly_list)
    
    display(plot([power_range1 power_range2], [times1 times2], 
        title = "Comparing methods with varying power",
        label = ["New" "Old"], 
        xlabel = "Power", 
        ylabel = "Multiplication time (ms)", 
        shape = :circle,
        legend = :topleft))

end

plot_pow_timings();