
using NLPModels
#using CUTEst
using OptimizationProblems
using NLPModelsJuMP
using SolverTools   # Pour avoir les utilitaires d'affichage log_header et log_row

using Stopping
using Logging


#include("params.jl")
γ = 1.0e-5
a = 0.0
b = Inf
scale = 1.0e-5

#nlp = CUTEstModel("BEALE")
nlp = MathOptNLPModel(beale())

dir = -grad(nlp, nlp.meta.x0) * scale
h = LineModel(nlp, nlp.meta.x0, dir);
dϕt0 = grad(h,a)

ϵ = max(abs(dϕt0),1.0)
α = -γ*ϵ
β =  γ*ϵ

τ₁ = 0.99
τ₀ = 0.001
g₀ = dϕt0


#logger = Logging.ConsoleLogger(stderr,Logging.Warn)
#logger = Logging.ConsoleLogger(stderr,Logging.Info)
logger = Logging.ConsoleLogger(stderr,Logging.Debug)
#logger = Logging.NullLogger()


stp = NLPStopping(h, OneDAtX(0.0,zeros(2)), tol_check=(a,b,c) -> ones(2))
stp.meta.max_iter = 100
stp.meta.atol = γ
stp.meta.rtol = γ


stp = Glob_U(h, a=0.0, b=Inf, α=α, β =β, pick_in = pick_inN2, stp=stp)
@show stp.current_state.x, stp.current_state.fx, stp.current_state.gx, status(stp)
@show nlp.counters
@show h.counters

@test true

