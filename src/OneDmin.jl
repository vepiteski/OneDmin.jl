module OneDmin

using Polynomials
using LinearAlgebra

include("Local/LocalIterations.jl")

using NLPModels
using SolverTools
using SolverCore
using Logging
using Stopping

include("LineSearch.jl")

include("bracket/Bracket.jl")

include("TR/TR1D.jl")


end
