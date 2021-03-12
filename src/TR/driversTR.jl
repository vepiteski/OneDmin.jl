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

function TR1D_N(ϕ; kwargs...)
    return TR1D(ϕ; M = Taylor2{Float64}(), kwargs...)
end

function TR1D_s(ϕ; kwargs...)
    return TR1D(ϕ; M = secant{Float64}(), kwargs...)
end

function TR1D_C2(ϕ; kwargs...)
    return TR1D(ϕ; M = TwoPoints3{Float64}(), kwargs...)
end

function TR1D_s2(ϕ; kwargs...)
    return TR1D(ϕ; M = secant3{Float64}(), kwargs...)
end

