export TR1D_N, TR1D_s, TR1D_s2, TR1D_C2


""" `TR1D_N(ϕ, a, b, stp, α, β, η₁, η₂, red. aug, M)`
   - ϕ is a linemodel
   - kwargs from `bracket`
      - minimization in the interval [a,b]
      - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<ϕ(a)
      - stp : stopping object; atol and rtol are bypassed and the conditions above apply


`TR1D_N` is a scalar trust region implementation using the Taylor quadratic model (Newton's method)
"""
function TR1D_N(ϕ; kwargs...)
    return TR1D(ϕ; M = Taylor2{Float64}(), kwargs...)
end

""" `TR1D_N(ϕ, a, b, stp, α, β, η₁, η₂, red. aug, M)`
   - ϕ is a linemodel
   - kwargs from `bracket`
      - minimization in the interval [a,b]
      - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<ϕ(a)
      - stp : stopping object; atol and rtol are bypassed and the conditions above apply


`TR1D_s` is a scalar trust region implementation using the two-points secant  quadratic model
"""
function TR1D_s(ϕ; kwargs...)
    return TR1D(ϕ; M = secant{Float64}(), kwargs...)
end


""" `TR1D_C2(ϕ, a, b, stp, α, β, η₁, η₂, red. aug, M)`
   - ϕ is a linemodel
   - kwargs from `bracket`
      - minimization in the interval [a,b]
      - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<ϕ(a)
      - stp : stopping object; atol and rtol are bypassed and the conditions above apply


`TR1D_C2` is a scalar trust region implementation using the two-points cubic model. Convergence not proved.
"""
function TR1D_C2(ϕ; kwargs...)
    return TR1D(ϕ; M = TwoPoints3{Float64}(), kwargs...)
end


""" `TR1D_C2(ϕ, a, b, stp, α, β, η₁, η₂, red. aug, M)`
   - ϕ is a linemodel
   - kwargs from `bracket`
      - minimization in the interval [a,b]
      - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<ϕ(a)
      - stp : stopping object; atol and rtol are bypassed and the conditions above apply


`TR1D_s2` is a scalar trust region implementation using the two-points cubic model approximated by 2 secant sub-iterations. Convergence not proved.
"""
function TR1D_s2(ϕ; kwargs...)
    return TR1D(ϕ; M = secant3{Float64}(), kwargs...)
end

