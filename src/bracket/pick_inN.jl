#export pick_inN2, pick_inN, pick_inN3

"""
pick\\_inN2 is the picking scheme which picks the Newton point whenever possible. 

"""
function pick_inN2(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
    γ = 0.9
    L = b - a
    if forward
        dN = - (dϕa / hess(ϕ, a))
        if (dN > 0) & (dN < γ*L)
            return a + dN             
        elseif b<Inf
            return (a + b)/2.0
        else
            return pick_inInf(a)
        end
    else
        dN =  - (dϕb / hess(ϕ, b))
        if (dN < 0) & (dN > -γ*L)
            return b + dN
        else
            return (a + b)/2.0
        end
    end
end


"""
pick\\_inN is the picking scheme which picks the Newton point whenever possible but 
         uses pick\\_inInf whenever b=∞

Pick\\_inN2 is more efficient.

"""
function pick_inN(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
    γ = 0.9
    if b<Inf
        L = b - a
        if forward
            dN = - (dϕa / hess(ϕ, a))
            if (dN > 0) & (dN < γ*L)
                return a + dN             
            else
                return (a + b)/2.0
            end
        else
            dN =  - (dϕb / hess(ϕ, b))
            if (dN < 0) & (dN > -γ*L)
                return b + dN
            else
                return (a + b)/2.0
            end
        end
    else
        return pick_inInf(a)
    end
end


# Tricky implementation to keep the predecessor iterate for inclusion in the secant formula
#
# 

let (tpred, ϕpred, dϕpred, d2ϕpred) = (Inf, Inf, Inf, Inf)
    γ = 0.9
    global function resetN3!(;t::Float64 = Inf, vϕ::Float64 = Inf, dϕ::Float64 = Inf, d2ϕ::Float64 = Inf)
        (tpred, ϕpred, dϕpred, d2ϕpred) = (t, vϕ, dϕ, d2ϕ)
    end

    # based on pick_inS2Opt: filters first the first Newton step and then tries the second.    

"""
pick\\_inN3 is the picking scheme which picks the Newton point whenever possible and then eventually
     performs a further pseudo-Newton iteration using a quintic Hermite interpolant

"""
    global function pick_inN3(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
        γ = 0.9
        L = b - a
        if forward
            d2ϕa = hess(ϕ, a)
            dN = - (dϕa / d2ϕa)
            t2 = a + dN
            # first the Newton's improvement
            if (dN > 0) & (dN < γ*L)
                #prepare to improve the Newton step
                α = coefHIN5([a;tpred],[ϕa dϕa d2ϕa; ϕpred dϕpred d2ϕpred])
                p = Polynomial(α)
                dp = Polynomials.derivative(p)
                d2p = Polynomials.derivative(dp)
                d2 = -(dp(dN) / d2p(dN))
                dN2 = dN + d2
                d3 = -(dp(dN2) / d2p(dN2))
                dN3 = dN2 + d3
                d4 = -(dp(dN3) / d2p(dN3))
                dN4 = dN3 + d4
               if (dN4 > 0) & (dN4 < γ*L)
                    t2 = a + dN4
                end     
            elseif b<Inf
                t2 = (a + b)/2.0
            else
                t2 = pick_inInf(a)
            end
            tpred, ϕpred, dϕpred, d2ϕpred = a, ϕa, dϕa, d2ϕa
            return t2
        else
            d2ϕb = hess(ϕ, b)
            dN =  - (dϕb / d2ϕb)
            t2 = b + dN
            # first the Newton's improvement
            if (dN < 0) & (dN > -γ*L)
                 #prepare to improve the Newton step
                α = coefHIN5([b;tpred],[ϕb dϕb d2ϕb; ϕpred dϕpred d2ϕpred])
                p = Polynomial(α)
                dp = Polynomials.derivative(p)
                d2p = Polynomials.derivative(dp)
                d2 = -(dp(dN) / d2p(dN))
                dN2 = dN + d2
                d3 = -(dp(dN2) / d2p(dN2))
                dN3 = dN2 + d3
                d4 = -(dp(dN3) / d2p(dN3))
                dN4 = dN3 + d4
                if (dN4 > 0) & (dN4 < γ*L)
                    t2 = b + dN4
                end     
            else
                t2 = (a + b)/2.0
            end
            tpred, ϕpred, dϕpred, d2ϕpred = b, ϕb, dϕb, d2ϕb
            return t2
        end
    end
    
    
end




"""
pick\\_inN3 is the picking scheme which picks the Newton point whenever possible and then eventually
     performs a further pseudo-Newton iteration using a quintic Hermite interpolant

"""
pick_inN3
