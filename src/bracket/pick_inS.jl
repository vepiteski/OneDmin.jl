export pick_ins2Opt, pick_ins, pick_ins2, pick_ins2N, pick_ins1

# Tricky implementation to keep the predecessor iterate for inclusion in the secant formula
#
# 

let (tpred, ϕpred, dϕpred) = (Inf, Inf, Inf)
    γ = 0.7
    γ2 = 0.999
    global function resetSec!(;t::Float64 = Inf, vϕ::Float64 = Inf, dϕ::Float64 = Inf)
        (tpred, ϕpred, dϕpred) = (t, vϕ, dϕ)
    end

    # test version that filters first the first secant step and then tries the second.
    # Should give the same iterates, slightly faster.
    
    """
        pick_ins2Opt is the scheme which picks the secant point whenever possible and improves it
        using a pseudo further secant iteration using an Hermite interpolant of ϕ using ϕand dϕ at 
        the points defining the secant.
    """
    global function pick_ins2Opt(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)        
        L = b - a
        if forward
            s = tpred - a
            y = dϕpred - dϕa
            dS = - (dϕa * s / y)
            t2 = a + dS
            # first check secant improvement
            if (dS > 0) & (dS < γ*L)
                # prepare to improve the secant step
                α = coefHIN3([a;tpred],[ϕa dϕa; ϕpred dϕpred])
                p = Polynomial(α)
                dp = Polynomials.derivative(p)
                dϕ2 = dp(dS)
                y2 = dϕ2 - dϕa
                d2 = - (dϕ2 * dS / y2)
                dS2 = dS + d2
                if (dS2 > 0) & (dS2 < γ*L)
                    t2 = a + dS2
                else
                end
            elseif b<Inf
                t2 = (a + b)/2.0
            else
                t2 = pick_inInf(a)
            end
            tpred, ϕpred, dϕpred = a, ϕa, dϕa
            return t2
        else
            s = tpred - b
            y = dϕpred - dϕb
            dS =  - (dϕb * s / y)
            t2 = b + dS
            # first check secant improvement
            if (dS < 0) & (dS > -γ*L)
                # prepare to improve the secant step
                α = coefHIN3([b;tpred],[ϕb dϕb; ϕpred dϕpred])
                p = Polynomial(α)
                dp = Polynomials.derivative(p)
                dϕ2 = dp(dS)
                y2 = dϕ2 - dϕb
                d2 = - (dϕ2 * dS / y2)
                dS2 = dS + d2
                if (dS2 < 0) & (dS2 > -γ*L)
                    t2 = b + dS2
                else 
                    t2 = b + dS
                end
            else
                t2 = (a + b)/2.0
            end
            tpred, ϕpred, dϕpred = b, ϕb, dϕb
            return t2
            
        end
    end
    

                              
    """
    pick\\_ins is the picking scheme which picks the secant point whenever possible. 
    """
    global function pick_ins(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
        if b < Inf
            L = b - a
            if forward
                s = tpred - a
                y = dϕpred - dϕa
                dS = - (dϕa * s / y)
                tpred, dϕpred = a, dϕa
                anext = a + dS
                if (dS > 0) & (dS < γ*L)
                    return anext
                else
                    return (a + b)/2.0
                end
            else
                s = tpred - b
                y = dϕpred - dϕb
                dS =  - (dϕb * s / y)
                tpred, dϕpred = b, dϕb
                bnext = b + dS
                if (dS < 0) & (dS > -γ*L)
                    return bnext
                else
                    return (a + b)/2.0
                end
            end
        else
            return pick_inInf(a)
        end
    end
    

    """
            pick\\_ins2 is the picking scheme which picks the secant point whenever possible and improves it
        using a pseudo further secant iteration using an Hermite interpolant of ϕ using ϕ and dϕ at 
        the points defining the secant.
    """
    global function pick_ins2(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
        L = b - a
        if forward
            s = tpred - a
            y = dϕpred - dϕa
            dS = - (dϕa * s / y)
            t2 = a + dS
            α = coefHIN3([a;tpred],[ϕa dϕa; ϕpred dϕpred])
            p = Polynomial(α)
            dp = Polynomials.derivative(p)
            dϕ2 = dp(dS)
            y2 = dϕ2 - dϕa
            d2 = - (dϕ2 * dS / y2)
            dS2 = dS + d2
            
            anext = a + dS2
            
            tpred, ϕpred, dϕpred = a, ϕa, dϕa
            if (dS2 > (1-γ2)*L) & (dS2 < γ*L)
                return anext
            elseif (dS > (1-γ2)*L) & (dS < γ*L)
                return a + dS
            elseif b<Inf
                return (a + b)/2.0
            else
                return pick_inInf(a)
            end
        else
            s = tpred - b
            y = dϕpred - dϕb
            dS =  - (dϕb * s / y)
            t2 = b + dS

            α = coefHIN3([b;tpred],[ϕb dϕb; ϕpred dϕpred])
            p = Polynomial(α)
            dp = Polynomials.derivative(p)
            dϕ2 = dp(dS)

            y2 = dϕ2 - dϕb
            d2 = - (dϕ2 * dS / y2)
            dS2 = dS + d2

            bnext = b + dS2

            tpred, ϕpred, dϕpred = b, ϕb, dϕb
            if (dS2 < (1-γ2)*L) & (dS2 > -γ*L)
                return bnext
            elseif (dS < (1-γ2)*L) & (dS > -γ*L)
                return b + dS
            else
                return (a + b)/2.0
            end
        end
    end


    
    """
        pick\\_ins2N is the picking scheme picking the secant point whenever possible and improves it
        using a pseudo further Newton iteration using an Hermite interpolant of ϕ using ϕ and dϕ at 
        the points defining the secant.
    """
    global function pick_ins2N(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
        L = b - a
        if forward
            s = tpred - a
            y = dϕpred - dϕa
            dS = - (dϕa * s / y)
            t2 = a + dS
            α = coefHIN3([a;tpred],[ϕa dϕa; ϕpred dϕpred])
            p = Polynomial(α)
            dp = Polynomials.derivative(p)
            d2p = Polynomials.derivative(dp)
            dϕ2 = dp(dS)
            d2ϕ2 = d2p(dS)
            d2 = - (dϕ2 / d2ϕ2)
            dS2 = dS + d2
            
            anext = a + dS2
            
            tpred, ϕpred, dϕpred = a, ϕa, dϕa
            if (dS2 > 0) & (dS2 < γ*L)
                return anext
            elseif (dS > 0) & (dS < γ*L)
                return a + dS
            elseif b<Inf
                return (a + b)/2.0
            else
                return pick_inInf(a)
            end
        else
            s = tpred - b
            y = dϕpred - dϕb
            dS =  - (dϕb * s / y)
            t2 = b + dS

            α = coefHIN3([b;tpred],[ϕb dϕb; ϕpred dϕpred])
            p = Polynomial(α)
            dp = Polynomials.derivative(p)
            d2p = Polynomials.derivative(dp)
            dϕ2 = dp(dS)

            d2ϕ2 = d2p(dS)
            d2 = - (dϕ2 / d2ϕ2)
            dS2 = dS + d2

            bnext = b + dS2

            tpred, ϕpred, dϕpred = b, ϕb, dϕb
            if (dS2 < 0) & (dS2 > -γ*L)
                return bnext
            elseif (dS < 0) & (dS > -γ*L)
                return b + dS
            else
                return (a + b)/2.0
            end
        end
    end

    """
        pick\\_ins2 is the picking scheme which picks the secant point whenever possible and improves it
        using a pseudo further secant iteration using an Hermite interpolant of ϕ using ϕ and dϕ at 
        the points defining the secant.
    """
    global function pick_ins1(a, b, ϕ, ϕa, ϕb, dϕa, dϕb, forward)
        L = b - a
        if forward
            s = tpred - a
            y = dϕpred - dϕa
            dS = - (dϕa * s / y)
            t = a + dS
            
            tpred, ϕpred, dϕpred = a, ϕa, dϕa
            if (dS > (1-γ2)*L) & (dS < γ*L)
                return t
            elseif b<Inf
                return (a + b)/2.0
            else
                return pick_inInf(a)
            end
        else
            s = tpred - b
            y = dϕpred - dϕb
            dS =  - (dϕb * s / y)
            t = b + dS
            tpred, ϕpred, dϕpred = b, ϕb, dϕb
            if (dS < (1-γ2)*L) & (dS > -γ*L)
                return t
            else
                return (a + b)/2.0
            end
        end
    end


end  #  let

# Repeat the doc, the actual behavior of docstrings is not as documented:(


    """
    pick\\_ins is the picking scheme which picks the secant point whenever possible. 
    """
pick_ins


    """
        pick\\_ins2N is the picking scheme picking the secant point whenever possible and improves it
        using a pseudo further Newton iteration using an Hermite interpolant of ϕ using ϕ and dϕ at 
        the points defining the secant.
    """
pick_ins2N


    """
        pick_ins2Opt is the scheme which picks the secant point whenever possible and improves it
        using a pseudo further secant iteration using an Hermite interpolant of ϕ using ϕand dϕ at 
        the points defining the secant.
    """
pick_ins2Opt


    """
        pick\\_ins1 is the picking scheme which picks the secant point whenever possible and improves it
        using a pseudo further secant iteration using an Hermite interpolant of ϕ using ϕ and dϕ at 
        the points defining the secant.
    """
pick_ins1


"""
pick\\_ins2 is the picking scheme which picks the secant point whenever possible and improves it
        using a pseudo further secant iteration using an Hermite interpolant of ϕ using ϕ and dϕ at 
        the points defining the secant.
"""
pick_ins2



