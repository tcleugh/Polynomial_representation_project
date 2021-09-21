



Random.seed!(0)

function do_timing()
    data = rand(Int,3*10^4)
    prediction_seconds =  predict_bubble_time(3*10^4)/10^3
    print("Predicted time: $prediction_seconds \n Actual time:")
    @time bubble_sort!(data);

    data = rand(Int,5*10^4)
    prediction_seconds =  predict_bubble_time(5*10^4)/10^3
    print("\nPredicted time: $prediction_seconds \n Actual time:")
    @time bubble_sort!(data);
end
do_timing();