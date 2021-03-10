using Plots
pyplot()   
using SolverBenchmark
using DataFrames

# We consider as solved the really solved (first_order) but also
# unbounded parameters (infeasible) or function values (unbounded)
# as well as exceptions (NaN in the evaluation of the model).
 

solved(df) = (df.status .== :first_order) .| (df.status .== :infeasible) .| (df.status .== :unbounded).| (df.status .== :exception)
 


costnames = ["time",
             "objective evals",
             "derivative evals",
             "second derivative evals",
             "obj + der + 2nd der"]

costs = [df -> .!solved(df) .* Inf .+ df.elapsed_time,
         df -> .!solved(df) .* Inf .+ df.neval_obj,
         df -> .!solved(df) .* Inf .+ df.neval_grad,
         df -> .!solved(df) .* Inf .+ df.neval_hess,
         df -> .!solved(df) .* Inf .+ df.neval_obj .+ df.neval_grad .+ df.neval_hprod]

#include("my_profile_tools.jl")

# profile on the sum of evaluations obj, grad and hess
#p = my_performance_profile(stats,costs[1], title = costnames[1])

# array of profiles
#p = my_profile_solvers(pruned_stats, costs, costnames)


#Plots.pdf(p, "test") 
