include("TR1DModels.jl")
include("TR1DModels3.jl")
include("driversTR.jl")

export TR1D

""" `TR1D(ϕ, a, b, stp, α, β, η₁, η₂, red. aug, M)`
   - ϕ is a linemodel
   - minimization in the interval [a,b], default [0,∞]. 
      - ϕ'(a) is assumed negative
   - stp : stopping object; atol and rtol are bypassed and the conditions below apply
   - α<0<β is the (possibly asymmetric) stopping tolerance, default β=-α=1e-8
      - solution to satisfy   α ≤ ϕ'(t\\*) ≤ β and ϕ(t\\*)<h(a)
   - η₁ and η₂ : trust region treshold for unsuccessesful and very successful iterations, default η₁=0.25, η₂ = 0.75
   - red, aug : trust region reduction and augmentation ratio, default red=0.5, aug=5
   - M a model in the sense of trust region, usually quadratic Taylor model.

`TR1D` is a scalar trust region implementation
"""
function TR1D(ϕ       :: OneDModel;
              a       :: T = 0.0, b :: T = Inf,
              stp     :: AbstractStopping =
              NLPStopping(ϕ, OneDAtX(a, zeros(2)),
                          tol_check=(a,b,c) -> ones(2)),
              #stp     :: AbstractStopping = LS_Stopping(ϕ, LSAtT(a)),
              α       :: T = -1e-8,  β  :: T = 1e-8,
              η₁      :: Float64 = 0.25,
              η₂      :: Float64 = 0.75,
              red     :: Float64 = 0.5,
              aug     :: Float64 = 5.0,
              M       :: TR1DModel = Taylor2(0.0,0.0,0.0) )  where T
    # Specialized TR for handling non-negativity constraint on t, i.e. a <= t
    # Trust region parameters
    Δp = 1.0  # >=0
    Δn = 0.0  # <=0
    
    t=a
    q(d) = m(M, d)
        
    ϕt = obj(ϕ, t)
    dϕt = grad(ϕ,t)
    dϕtd = NaN
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
    stp.meta.tol_check     = (a,b,c) -> [  ϕt; β]
    
    M = initModel(M, ϕ, t, ϕt, dϕt)
    @info log_header([:iter, :f, :dual, :step, :slope, :success], [Int, T, T, T, T, String],
                     hdr_override=Dict(:f=>"h(t)", :dual=>"h'(t)", :step=>"Δn", :slope=>"Δp"))

    #OK = update_and_start!(stp, x = a, fx=ϕt, gx = dϕt)
    ϕ1  = obj(ϕ, T(1))
    dϕ1 = grad(ϕ, T(1))
    OK = update_and_start!(stp, x = T(1), fx=ϕ1, gx = dϕ1) # try one for line search

    while !OK
        dN = solveModel(M, t)
        #if (q(Δp)<q(Δn)) #| (Δn==0.0)
        if (Δm(M, Δp, Δn) < 0.0) #| (Δn==0.0)
            d=Δp
            @debug "dp"
        else
            d=Δn
            @debug "dn"
        end
        
        #if  (dN < Δp) & (dN > Δn) & (q(d)>q(dN))
        if  (dN < Δp) & (dN > Δn) & (Δm(M, d, dN) > 0.0)
#        if  (dN < Δp) & (dN > Δn) & ( (dϕt+0.5*H*(d+dN))*(d-dN) >0.0)
#        if  (dN < Δp) & (dN > Δn) & ( Δm(M, d-dN) >0.0)
            d=dN 
            #@show "dN = ", d
        end
        
        ϕtestTR = obj(ϕ,t+d)
        pred = Δm(M, d)
        if pred > 0   @show pred  end
        @debug "d = ", d
        @debug "H = ", M.H
        @debug "dϕt = ", M.g
        @debug "pred = ", pred

        if (pred > 0) #|| (Δp < 1.0)
            @show pred
            @show M
            @show d, dN, Δp, Δn
            @show Δm(M, dN)
            @show q(d), q(dN), q(d)-q(dN)

            disc = M.H^2 - 2.0*M.C*M.g
            dC = disc < 0 ? -sign(M.C)*Inf : (-M.H + sqrt(disc)) / M.C
            @show dC
            
            tg = t+2.0*Δn:0.01*Δp:t+1.0*Δp
            n, = size(tg)
            fg=zeros(n)
            qg = zeros(n)
            dqg = zeros(n)
            for i = 1:n
                fg[i] = obj(ϕ,tg[i])
                #qg[i] = q(tg[i]-t)
                qg[i] = m(M, tg[i]-t)
                dqg[i] = dm(M, tg[i]-t)
            end
            @show "t = ", t, "d = ", d
            @show "Δp = ", Δp, "Δn = ", Δn
            plt = plot(tg,[fg, qg, dqg], label = ["f"  "p"  "dp"])
            return plt
            @assert pred<=0 
        end


        goodgrad = false
        #ared = ϕtestTR-ϕt # inclure ici le truc numérique de H & Z
        if abs(d) < 1e-8  dϕtd = grad(ϕ,t+d); ared = (dϕtd + dϕt)*d/2.0; goodgrad = true
        else  ared = ϕtestTR-ϕt
        end
        
        ratio = ared / pred 
        
        if ratio < η₁ # Unsuccessful
            @debug ared,  pred,  t,  d, dϕt 
            Δp = red*Δp
            Δn = red*Δn
            @info log_row(Any[stp.meta.nb_of_stop, ϕt, dϕt, Δn, Δp, "U"])
        else             # Successful
            t = t + d
            if goodgrad
                dϕt = dϕtd
            else
                dϕt = grad(ϕ,t)
            end
            
            ϕt = ϕtestTR

            M = updateModel!(M, ϕ, t, ϕt, dϕt)
            
            if ratio > η₂ 
                Δp = aug * Δp
                Δn = max(a, aug * Δn) - t
            else
                Δn = max(a, Δn) - t
            end
            OK = update_and_stop!(stp, x=t, fx=ϕt, gx=dϕt)
            @info log_row(Any[stp.meta.nb_of_stop, ϕt, dϕt, Δn, Δp,  "S"])
        end;
    end;
        
    return stp
    
end

