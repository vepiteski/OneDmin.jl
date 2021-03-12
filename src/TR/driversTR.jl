export TR1D_N, TR1D_s, TR1D_s2, TR1D_C2


""" `TR1D_N(ϕ, a, b, stp, α, β, η₁, η₂, red. aug, M)`
   - ϕ is a linemodel
   - minimization in the interval [a,b]
      - ϕ'(a) is assumed negative
   - α<0<β is the (possibly asymmetric) stopping tolerance, default β=-α=1e-8
      - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<h(a)
   - stp : stopping object; atol and rtol are bypassed and the conditions above apply

`TR1D_N` is a scalar trust region implementation using the Taylor quadratic model (Newton's method)
"""

function TR1D_N(ϕ, a, b, α, β, stp; kwargs...)
    return TR1D(ϕ, a=a, b=b, α=α, β=β, stp=stp; M = Taylor2{Float64}(), kwargs...)
end

function TR1D_s(h, a, b, α, β, stp; kwargs...)
    return TR1D(h, a=0.0, b=Inf, α=α, β=β, stp=stp; M = secant{Float64}(), kwargs...)
end

function TR1D_C2(h, a, b, α, β, stp; kwargs...)
    return TR1D(h, a=0.0, b=Inf, α=α, β=β, stp=stp; M = TwoPoints3{Float64}(), kwargs...)
end

function TR1D_s2(h, a, b, α, β, stp; kwargs...)
    return TR1D(h, a=0.0, b=Inf, α=α, β=β, stp=stp; M = secant3{Float64}(), kwargs...)
end

