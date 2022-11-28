# How to call?
# in a line search algo, where x is the current point, f current objective, g current grad
# and nlp the model

# A LSModel implements the ϕ function. An approximate minimizer of ϕ such that ϕ(t)<ϕ(0)
# and such that εa ≤ ϕ'(t) ≤ ɛb  satisfies the (strong) Wolfe conditions
# the ϕ function needs to know τ₀ and should store f₀ and g₀ to avoid recomputations

#export LSModel, OneDModel, obj, grad, grad!, hess, rebase!, linesearch, prepare_LS 
export LSModel, OneDModel, linesearch, prepare_LS

""" A LSModel is a specialized LineModel keeping track of the objective and gradient
    of the original NLPModel from which it is built. Moreover, it adjust the line model
    to the Armijo parameter τ₀: ϕ(t) = f(x₀ + t*d) - f(x₀) - t*τ₀*∇f(x₀)*d
    """
mutable struct LSModel{T, S} <: AbstractNLPModel{T, S}
    meta :: NLPModelMeta{T, S}
    counters :: Counters
    nlp :: AbstractNLPModel{T, S}
    x₀:: S
    d :: S
    τ₀:: T
    f₀:: T   # f(x₀)
    g₀:: T   # ∇f(x₀)⋅d
    x :: S   # current point in the nlp
    f :: T   # objective of the nlp f(x)
    ∇f:: S   # to avoid recomputing it outside the line search 
end


OneDModel = Union{LineModel, LSModel}



function LSModel(nlp::AbstractNLPModel{T, S}, x₀::S, d::S, τ::T) where {T, S}
    f₀ = obj(nlp, x₀)
    g₀ = grad(nlp, x₀)⋅d
    ∇f = similar(x₀)   # instantiate the room to memorize the last gradient
    x = copy(x₀)
    τ₀ = τ
    meta = NLPModelMeta{T, S}(1, x0=zeros(T, 1), name="LSModel to $(nlp.meta.name))")

    return LSModel(meta, Counters(), nlp, x₀, d, τ₀, f₀, g₀, x, f₀, ∇f)
end

function LSModel(nlp::AbstractNLPModel{T, S}, x₀::S, d::S, τ::T, f₀::T, g₀::T) where {T, S}
    ∇f = similar(x₀)   # instantiate the room to memorize the last gradient
    x = copy(x₀)
    τ₀ = τ
    meta = NLPModelMeta(1, x0=zeros(T, 1), name="LSModel to $(nlp.meta.name))")

    return LSModel(meta, Counters(), nlp, x₀, d, τ₀, f₀, g₀, x, f₀, ∇f)
end



import NLPModels: obj, grad, grad!, hess

function obj(ϕ::LSModel, t::AbstractFloat)
    NLPModels.increment!(ϕ, :neval_obj)
    ϕ.x = ϕ.x₀ + t*ϕ.d
    ϕ.f = obj(ϕ.nlp, ϕ.x)
    f =  ϕ.f - ϕ.f₀ - t*ϕ.τ₀*ϕ.g₀
    return f
end


function grad(ϕ::LSModel, t::AbstractFloat)  # stocks the gradient of the nlp in ϕ.∇f
    NLPModels.increment!(ϕ, :neval_grad)
    ϕ.x = ϕ.x₀ + t*ϕ.d
    g =  grad!(ϕ.nlp, ϕ.x, ϕ.∇f)⋅ϕ.d - ϕ.τ₀*ϕ.g₀
    return g
end


function grad!(ϕ :: LSModel, t :: AbstractFloat, gv :: AbstractVector)
    NLPModels.increment!(ϕ, :neval_grad)
    ϕ.x = ϕ.x₀ + t*ϕ.d
    iszero(t) ? (1-ϕ.τ₀)*g₀ : grad!(ϕ.nlp, ϕ.x, gv)⋅ϕ.d - ϕ.τ₀*ϕ.g₀
    return g

end


function hess(ϕ::LSModel, t::AbstractFloat)
    NLPModels.increment!(ϕ, :neval_hess)
    ϕ.x = ϕ.x₀ + t*ϕ.d
    H = hprod(ϕ.nlp,  ϕ.x, ϕ.d)⋅ϕ.d
    return H
end



"""`rebase!(ϕ, x, d, τ₀, f₀, g₀)`
Change the values of x₀, d, τ₀, of the LSModel ϕ, but retains the counters.
    If f₀ and/or g₀ are not provided, they are computed to set ϕ.f₀ and ϕ.g₀.
"""
function rebase!(ϕ :: LSModel, x₀ :: AbstractVector, d :: AbstractVector, τ₀ :: AbstractFloat,
                 f₀ :: AbstractFloat = nothing, g₀ :: AbstractFloat = nothing)
    ϕ.x₀, ϕ.d, ϕ.τ₀ = x₀, d, τ₀
    if f₀ == nothing f₀ = obj(ϕ.nlp,x₀); end
    ϕ.f₀ = f₀
    if g₀ == nothing g₀ = dot(grad(ϕ.nlp, x₀),d); end
    ϕ.g₀ = g₀
    
    return ϕ
end

# With that given, the setup for a linesearch within an iteration is to call linesearch within
# the loop, and prepare_LS before the loop.

function linesearch(ϕ    :: LSModel,
                    ϕstp :: AbstractStopping,  # be more specific for stp
                    x₀   :: Vector{T},
                    d    :: Vector{T},
                    f₀   :: T,
                    ∇f₀  :: Vector{T},
                    τ₀   :: Real,
                    τ₁   :: Real;
                    strongWolfe :: Bool = false,
                    logger :: AbstractLogger = NullLogger(),
                    algo :: Function = bracket,
                    kwargs...) where T<:AbstractFloat

    # rebase the LSModel 
    g₀ = dot(∇f₀,d)
    rebase!(ϕ, x₀, d, τ₀, f₀, g₀)
    
    # convert the Armijo and Wolfe criteria to an asymetric interval [α,β]
    α = T((τ₁-τ₀)*g₀)
    β = T(Inf) 
    if strongWolfe  β = -(τ₁+τ₀)*g₀ end

    # reuse the stopping
    reinit!(ϕstp, rlist = false)
    ϕstp.pb = ϕ
    # define the optimality_check function using α and β
    # ϕstp.meta.optimality_check =  (p,s) -> optim_check_1D(p,s,α,β)

    # optimize
    @debug "dans LS", α, β

    Logging.with_logger(logger) do
        #ϕstp = Glob_U(ϕ, a=0.0, b=Inf, α=α, β=β, stp = ϕstp, pick_in = pick_ins2, best=true)
        ϕstp = algo(ϕ, a=T(0.0), b=T(Inf), α=α, β=β, stp = ϕstp; kwargs...)
    end

    if !ϕstp.meta.optimal
        @warn "LineSearch failure"
        @info status(ϕstp,list=true)
    end
    # unpack the results
    t = ϕstp.current_state.x
    xt = ϕ.x
    ft = ϕ.f   # must rely on the solver (Glob_U) so that the last evaluation was
    gt = ϕ.∇f  # at the optimal step, so that the stored value for x, f and ∇f are valid

    return t, xt, ft, gt
end


function prepare_LS(stp :: AbstractStopping,
                    x₀  :: Vector{T},
                    d   :: Vector{T},
                    τ₀  :: Real,
                    f₀  :: Real,
                    g₀  :: Vector{T}) where T<:AbstractFloat

    # extract the nlp
    nlp = stp.pb
    
    # construct the line search model, which will be rebased at each iteration'current data
    ϕ = LSModel(nlp, x₀, d, τ₀, f₀, g₀⋅d)

    # instantiate stp, which will be adjusted at each iteration's current data
    #ϕstp = LS_Stopping(ϕ, (p,s) -> optim_check_LS(p,s), LSAtT(0.0),  main_stp = stp)
    ϕstp = NLPStopping(ϕ, OneDAtX(T(0.0),zeros(T,2)),
                       main_stp = stp,
                       tol_check=(a,b,c) -> ones(T,2))
    #ϕstp.stop_remote = cheap_stop_remote_control(resources_check=false)

    #
    #Optimality_check will be setup in the LineSearch call
    #

    ϕstp.meta.max_iter = -Int(log2(eps(T)))   #  52 Float64 precision 2^(-52)
    
    ϕstp.meta.atol = 0.0   # to rely only on the Armijo-Wolfe conditions
    ϕstp.meta.rtol = 0.0   # otherwise, tolerance may prohibe convergence
    ϕstp.meta.unbounded_threshold = 1e100
    return ϕ, ϕstp
end
