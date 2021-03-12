pickN = Function[]


""" `bracket_N(ϕ, a, b, α, β, stp;  kwargs...)`
   - ϕ is a linemodel
   - minimization in the interval [a,b]
   - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<ϕ(a)
   - stp : stopping object; atol and rtol are bypassed and the conditions above apply

driver to call bracket with pick_inN, the plain Newton pick in scheme 
"""
function bracket_N(ϕ, a, b, α, β, stp;  kwargs...)
    reset!(ϕ)
    return bracket(ϕ, a=a, b=b, α=α, β=β, stp = stp,  pick_in = pick_inN; kwargs...)
end
push!(pickN, bracket_N)


""" `bracket_N2(ϕ, a, b, α, β, stp;  kwargs...)`
   - ϕ is a linemodel
   - minimization in the interval [a,b]
   - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<ϕ(a)
   - stp : stopping object; atol and rtol are bypassed and the conditions above apply

driver to call bracket with pick_inN2, the variant Newton pick in scheme 
"""
function bracket_N2(ϕ, a, b, α, β, stp;  kwargs...)
    reset!(ϕ)
    return bracket(ϕ, a=a, b=b, α=α, β=β, stp = stp,  pick_in = pick_inN2; kwargs...)
end
push!(pickN, bracket_N2)


""" bracket_N3(h, a, b, α, β, stp;  kwargs...)
   - ϕ is a linemodel
   - minimization in the interval [a,b]
   - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<ϕ(a)
   - stp : stopping object; atol and rtol are bypassed and the conditions above apply

driver to call bracket with pick_inN3, the accelerated Newton pick in scheme 
"""
function bracket_N3(ϕ, a, b, α, β, stp;  kwargs...)
    reset!(ϕ)
    resetN3!()
    return bracket(ϕ, a=a, b=b, α=α, β=β, stp = stp,  pick_in = pick_inN3; kwargs...)
end
push!(pickN, bracket_N3)

picks = Function[]



""" driver to call bracket with pick_ins, the plain secant pick in scheme """
function bracket_s(h, a, b, α, β, stp;  kwargs...)
    reset!(h)
    resetSec!()
    return bracket(h, a=a, b=b, α=α, β=β, stp = stp,  pick_in = pick_ins; kwargs...)
end
push!(picks, bracket_s)


""" driver to call bracket with pick_ins2, the accelerated secant pick in scheme 
using 2 secant iterations on the cubic interpolant
"""
function bracket_s2(h, a, b, α, β, stp;  kwargs...)
    reset!(h)
    resetSec!()
    return bracket(h, a=a, b=b, α=α, β=β, stp = stp,  pick_in = pick_ins2; kwargs...)
end
push!(picks, bracket_s2)


""" driver to call bracket with pick_ins2N, the accelerated secant pick in scheme 
using 2 Newton iterations on the cubic interpolant
"""
function bracket_s2N(h, a, b, α, β, stp;  kwargs...)
    reset!(h)
    resetSec!()
    return bracket(h, a=a, b=b, α=α, β=β, stp = stp,  pick_in = pick_ins2N; kwargs...)
end
push!(picks, bracket_s2N)


""" driver to call bracket with pick_ins2Opt, the accelerated secant pick in scheme 
using 2 Newton iterations on the cubic interpolant (variant)
"""
function bracket_s2Opt(h, a, b, α, β, stp;  kwargs...)
    reset!(h)
    resetSec!()
    return bracket(h, a=a, b=b, α=α, β=β, stp = stp,  pick_in = pick_ins2Opt; kwargs...)
end
push!(picks, bracket_s2Opt)


""" driver to call bracket with pick_insB, the basic bissection pick in scheme 
"""

function bracket_B(h, a, b, α, β, stp;  kwargs...)
    reset!(h)
    return bracket(h, a=a, b=b, α=α, β=β, stp = stp,  pick_in = pick_inB; kwargs...)
end


export pickN, picks

for n in nameof.(pickN ∪ picks )
    @eval export $n
end

export bracket_B
