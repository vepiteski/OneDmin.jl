include("pick_inN.jl")
include("pick_inS.jl")
include("driversBracket.jl")


# Version inspired by JCG algo 6.6 (Fletcher-Lemaréchal) which avoids computing
# derivatives when not needed, e.g. if pick_in is the basic mid point.

export bracket


"""  `bracket(ϕ, a, b, stp, α, β, pick_in, best)`
   - ϕ is a linemodel
   - minimization in the interval [a,b], default [0,∞]. 
      - [a,b] is an initial interval quaranteed to contain a local minimum; if not satisfied, the function exits. In particular, ϕ'(a) must be negative.
   - stp : stopping object; atol and rtol are bypassed and the conditions below apply
   - α<0<β is the (possibly asymmetric) stopping tolerance, default β=-α=1e-8
   - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<h(a)
   - pick_in : function to pick a point in an interval; 
      -  several variants are proposed, default enhanced secant 
   - best option to bracket the best candidate or any minimizer such that h(t*)<h(a), default = true

`bracket` is a composite implementation of 
   - a bracketing scheme 
   - an interval reduction scheme using a pick_in function to select a point in a given interval
   - enhanced pick_in variants involving a polynomial interpolation (Newton like)
"""
function bracket(ϕ       :: OneDModel;   
                 a       :: T = Float64(0.0) ,   b  :: T = Inf,
                 #stp     :: AbstractStopping = NLPStopping(ϕ, OneDAtX(a)),
                 stp     :: AbstractStopping =
                            NLPStopping(ϕ, OneDAtX(a, zeros(T,2))),
                                        #tol_check=(a,b,c) -> ones(T,2)),
                 α       :: T = -1e-8,  β  :: T = 1e-8,
                 pick_in :: Function = pick_ins2,  # defaults to secant's enhancement
                 best    :: Bool = true,           # aim for the lowest minimizer
                 #kwargs...
                 ) where {T <: AbstractFloat}


    ϕa, dϕa = NaN,  NaN
    if isa(ϕ, LSModel) # then ϕ.f₀ and ϕ.g₀ are initialized
        @assert a==0
        ϕa, dϕa = T(0), (1-ϕ.τ₀)*ϕ.g₀
    else  # require a function and derivative computation 
        ϕa, dϕa = obj(ϕ, a), grad(ϕ, a)
    end

    tref = a
    ϕref = ϕa
    dϕref = dϕa

    
    ϕ₀ = ϕa
    ####################
    # à revoir, mais ça semble OK, on tente de trouver un meilleur minimiseur dans l'algo,
    # mais le critère d'arrêt ne teste que la descente par rapport à t=a.
    # on stoppe donc lorsque ϕ(t)<ϕ(a)   et   α ≤ dϕ(t) ≤  β
    #
    # Ce serait plus propre que ϕa=0 toujours...autrement, avec les erreurs d'arrondi,
    # le truc de tricher la valeur de ϕ courante ne marche possiblement pas... voir(*) plus bas
    ####################
    stp.meta.optimality_check = (p::OneDModel, s::AbstractState) -> [s.fx; s.gx] #
    stp.meta.atol = T(0)
    stp.meta.rtol = T(0)
    stp.meta.tol_check_neg = (a,b,c) -> [-Inf; α]
    stp.meta.tol_check     = (a,b,c) -> [  ϕ₀; β]


    
    ϕb, dϕb = NaN, NaN;    if isfinite(b)  ϕb, dϕb = obj(ϕ, b), grad(ϕ, b)  end

    tₘ, ϕₘ, dϕₘ = a, ϕa, dϕa # best objective so far

    OK, status = Check_coherent(a, b, ϕa, dϕa, ϕb, dϕb, α, β)
    if !OK
        @warn string("in one D minimization: ",status)
        @error "in one D minimization: non coherent"
        stp.meta.suboptimal = true
        return stp
    end

    OK = update_and_start!(stp, x = a, fx=ϕa, gx = dϕa)

    @info log_header([:iter, :f, :dual, :step, :slope], [Int, T, T, T, T],
                     hdr_override=Dict(:f=>"h(t)", :dual=>"h'(t)", :step=>"a", :slope=>"b"))

    @debug  "α = ", α, "   β = ", β
    dϕt = Inf
    forward = true
    initInf!()

    while !OK 
        t = pick_in(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
        # avoid computing derivatives as long as descent is not obtained
        ϕt = obj(ϕ, t); dϕt = Inf # grad(ϕ, t)

        # for as cheap a linesearch as possible, avoid to compare with ϕₘ, keep on
        #                                                 comparing to 0 (Armijo condition)
        if best          # rather compare with ϕₘ to aim the lowest minimizer possible.
            tref = tₘ
            ϕref = ϕₘ
            dϕref = dϕₘ
        end  
        
        descent = ϕt < ϕref
        if ~descent # double check for eventual roundoff error dominating the descent computation
            if (abs(ϕref - ϕt) < 1e-7) # heuristic, apply at roughly half the machine precision
                dϕt = grad(ϕ, t)       # not sure that the grad is already available
                val_descent = (dϕref+dϕt)*(t - tref)/2  
                    descent = val_descent < 0
            end
            # if descent, cheat and set the value of ϕt to its "descent" value       (*)
            if descent ϕt = ϕref + val_descent end
            @debug "best  ", ϕ₀, ϕₘ, ϕt, tₘ, t, dϕₘ, dϕt, descent
        end
        
        if descent
            # t and ϕt becomes the best guess so far
            tₘ, ϕₘ = t, ϕt
            dϕt = grad(ϕ, t)
            @debug " New grad dϕt = ", dϕt, "  t = ", t, "α = ", α, "   β = ", β
            @debug "           ϕt = ", ϕt, "  ϕ₀ = ", ϕ₀, " ϕₘ = ", ϕₘ
            dϕₘ = dϕt
            if dϕt > β
                b, ϕb, dϕb = t, ϕt, dϕt
                forward = false
            elseif dϕt < α
                a, ϕa, dϕa = t, ϕt, dϕt
                forward = true
            else  # should stop 
                @debug "   stop  dϕt =  $dϕt "
            end
        elseif forward
            b, ϕb = t, ϕt
        else  #backward
            a, ϕa = t, ϕt
        end
        @debug " avant stop   ϕt = $ϕt    ϕₘ = $ϕₘ, dϕₘ = $dϕₘ , dϕt = $dϕt "
        OK = update_and_stop!(stp, x=tₘ, fx=ϕₘ, gx=dϕₘ)
        @info log_row(Any[stp.meta.nb_of_stop, ϕt, dϕt, a, b])

    end
    if ϕₘ > ϕ₀ @warn "Increase in line minimization" end
    return stp
end



# Tricky implementation to have increasing jumps when b is infinite
let f_mult_inc = 5.0, inc = 1.0
    global function setInf!(v)
        inc = v
    end
    
    global function initInf!(;val :: Float64 = 1.0)
        setInf!(val)
        resetN3!() # resets the history for quintic Newton interpolation
        resetSec!()# resets the history for secant variants
    end

    global function getInf()
        return f_mult_inc, inc
    end
end

    
"""
    pick_inInf is the common increase function when b=Inf
"""
function pick_inInf(a)
    f_mult_inc, inc = getInf()
    anew = a + inc
    inc *= f_mult_inc
    setInf!(inc)
    return anew
end



"""
pick_inB is the basic picking scheme which picks the middle point in [a, b]. 

"""
function pick_inB(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
    
    if b < Inf
        return (a + b)/2.0
    else
        return pick_inInf(a)
    end
end


""" Check that the given interval [a,b] is guaranteed to contain a solution 
    such that α ≤ ϕ' ≤ β"""
function Check_coherent(a, b, ϕa, dϕa, ϕb, dϕb, α, β)

    OK = true
    status = :OK
    #  Check that the given interval [a,b] is guaranteed to contain a solution
    if b <= a
        OK = false
        status = :BadInterval
    elseif (dϕa >= 0.0)
        OK = false
        status = :NonDecrease
        #elseif !implies(b<Inf, (dϕb > β) || (ϕb >= ϕa))
        #    OK = false
        #    status = :NoSolutionGaranty
    elseif b<Inf
        if ϕb<ϕa
            if dϕb < β
                if dϕb > α
                    OK = true
                    status = :solved  #  b satisfies   (ϕb<ϕa)  &&  (α < dϕb < β)
                elseif dϕb < 0
                    OK = false
                    status = :NoSolutionWaranty
                end
            end
        end
    end

    return OK, status
end
