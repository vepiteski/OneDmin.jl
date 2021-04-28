using Pkg
Pkg.activate(".")



using Logging

using CUTEst
using NLPModels
using SolverBenchmark

using Printf
using DataFrames
using SolverTools
using Stopping

using OneDmin

probnames = sort(CUTEst.select(min_var=1,
                               max_var=1500,
                               max_con=0,
                               only_free_var=true))

@show length(probnames)
# Comment for the true test
 probnames = probnames[10:15]  # for testing quickly on a subset
# add a dummy copy of the first problem
probnames = insert!(probnames, 1, probnames[1])



problems = (CUTEstModel(probname) for probname in probnames)

γ = 1.0e-5
a = 0.0
b = Inf

scale = 1.0e-5

for (fn, fnsimple) ∈ [(:bracket_N, :N),(:bracket_s, :s),(:bracket_N3, :N3),(:bracket_s2Opt, :s2) ]

    @eval begin
        
        function $fnsimple(nlp)

            # Prepare a oneD model from a NLP model by using the line
            # model in the direction of the scaled gradient at the x0 point
            # of the NLP
            
            dir = -grad(nlp, nlp.meta.x0) * scale
            h = LineModel(nlp, nlp.meta.x0, dir);
            dϕt0 = grad(h,a)
            @assert dϕt0 < 0.0
            ϵ = max(abs(dϕt0),1.0)
            α = -γ*ϵ
            β =  γ*ϵ

            # setup the Stopping object
            stp = NLPStopping(h, OneDAtX(0.0,zeros(2)), tol_check=(a,b,c) -> ones(2))
            
            stp.meta.max_iter = 100

            reset!(h)
            reset!(nlp)
            reinit!(stp)

            OneDmin.initInf!()

            
            # Actual timing and execution
            t = @timed    stp = eval($fn)(h, a=a, b=b, α=α, β=β, stp=stp)
            
            # convert the information in Stopping to the "standard" diagnostic in JSO            
            st = stp.meta.optimal
            fb = !stp.meta.unbounded_pb
            tb = !stp.meta.unbounded
            ferr = stp.meta.domainerror
            tired = stp.meta.tired
            iter = stp.meta.nb_of_stop

            tsol = stp.current_state.x
            fh = stp.current_state.fx
            dh = stp.current_state.gx

            status = :unknown
            st &&    (status = :first_order)
            (!fb) && (status = :unbounded)
            (!tb) && (status = :infeasible)
            (ferr)&& (status = :exception)
            tired && (status = :max_iter)  # should be more specific
            
            
            return GenericExecutionStats(status, h, 
                                         solution = [tsol],
                                         iter = iter,
                                         primal_feas = tsol,  #  cheat to display the scalar solution
                                         dual_feas = abs(dh),
                                         objective = fh,
                                         elapsed_time = t[2],
                                         )
        end
        
    end
end


using OrderedCollections

solvers = OrderedDict{Symbol,Function}(
    :Sec => s,
    :Nwt => N,
    :Sec2 => s2,
    :Nwt3 => N3
)

# use versions accepting OrderedDict 
include("my_profile_tools.jl")

stats = my_bmark_solvers(solvers, problems,
                      colstats = [:name, :nvar, :status, :neval_obj, :neval_grad, :objective, :dual_feas, :primal_feas])

# remove the dummy copy of the first problem for statistics
for s in solvers
    stats[s[1]] = stats[s[1]][2:end,:]
end

include("PPpdfPaper.jl")

p = my_performance_profile(stats,costs[1], title = costnames[1])
Plots.pdf(p, "Figure2")
