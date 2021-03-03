module OneDmin

using Polynomials

include("LocalIterations.jl")

using NLPModels
using SolverTools
using Logging
using Stopping

include("LineSearch.jl")
include("pick_inN.jl")
include("pick_inS.jl")
include("Bracket.jl")



end
