include("TR1DModels.jl")
include("TR1DModels3.jl")
include("driversTR.jl")

export TR1D

"""
TR_U is a scalar trust region implementation allowing variants involving a polynomial 
interpolation (Newton like), in a Unified way

ϕ is the objective, a the starting point and it is assumed that ϕ'(a)<0 since 
only the half line t>=a is considered

A point t is searched such that ϕ(t) < ϕ(a)  and  α ≤ ϕ'(t) ≤ β
"""
function TR1D(ϕ       :: OneDModel;
              a       :: T = 0.0, b :: T = Inf,
              α       :: T = -1e-8,  β  :: T = 1e-8,
              stp     :: AbstractStopping = LS_Stopping(ϕ, LSAtT(a)),
              η₁      :: Float64 = 0.25,
              η₂      :: Float64 = 0.75,
              red     :: Float64 = 0.5,
              aug     :: Float64 = 5.0,
              verbose :: Bool = false,
              M       :: TR1DModel = secant(0.0,0.0,0.0) )  where T
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

function TR_UM_N(h, a, b, α, β, stp; kwargs...)
    return TR_UM(h, a=0.0, b=Inf, α=α, β=β, stp=stp; kwargs...)
end

function TR_UM_s(h, a, b, α, β, stp; kwargs...)
    return TR_UM(h, a=0.0, b=Inf, α=α, β=β, stp=stp; M = secant{Float64}(), kwargs...)
end

function TR_UM_C2(h, a, b, α, β, stp; kwargs...)
    return TR_UM(h, a=0.0, b=Inf, α=α, β=β, stp=stp; M = TwoPoints3{Float64}(), kwargs...)
end

function TR_UM_s2(h, a, b, α, β, stp; kwargs...)
    return TR_UM(h, a=0.0, b=Inf, α=α, β=β, stp=stp; M = secant3{Float64}(), kwargs...)
end


function TR_N(nlp, a, γ;  kwargs...)
    dir = -grad(nlp, nlp.meta.x0) * scale
    h = LineModel(nlp, nlp.meta.x0, dir);

    dϕt0 = grad(h,a)
    ϵ = max(abs(dϕt0),1.0)
    α = -γ*ϵ
    β =  γ*ϵ
    reset!(h)
    reset!(nlp)
    #@show a, b, α, β
    return TR_U(h, a, α, β, initH = initH_N, updateH = updateH_N ; kwargs...)
end


function TR_s(nlp, a, γ;  kwargs...)
    dir = -grad(nlp, nlp.meta.x0) * scale
    h = LineModel(nlp, nlp.meta.x0, dir);

    dϕt0 = grad(h,a)
    ϵ = max(abs(dϕt0),1.0)
    α = -γ*ϵ
    β =  γ*ϵ
    reset!(h)
    reset!(nlp)
    #@show a, b, α, β
    return TR_U(h, a, α, β, initH = initH_s!, updateH = updateH_s! ; kwargs...)
end


function TR_sec(nlp, a, b, γ, kwargs...)
    return TR_s(nlp, a, γ, kwargs...)
end

function TR_Nwt(nlp, a, b, γ; kwargs...)
    return TR_N(nlp, a, γ, kwargs...)
end


function TR_N_L(h, a, stp; verbose = false,  kwargs...)
    #@show "NEW TR_N"
    reset!(h)
    return TR_U(h, a, stp = stp,  M = Taylor2{Float64}(), verbose = verbose  ; kwargs...)
end

function TR_s_L(h, a, stp; verbose = false,  kwargs...)
    #@show "NEW TR_s"
    reset!(h)
    return TR_U(h, a, stp = stp,  M = secant{Float64}(), verbose = verbose  ; kwargs...)
end

function TR_cub2_L(h, a, stp; verbose = false,  kwargs...)
    #@show "NEW TR_cub2"
    reset!(h)
    return TR_U(h, a, stp = stp,  M = TwoPoints3{Float64}(), verbose = verbose ; kwargs...)
end

function TR_seccub2_L(h, a, stp; verbose = false,  kwargs...)
    #@show "NEW TR_seccub2"
    reset!(h)
    return TR_U(h, a, stp = stp,  M = secant3{Float64}(), verbose = verbose ; kwargs...)
end


function TR_sec_L(h, a, b, α, β, stp; verbose = false, kwargs...)
    return TR_s_L(h, a, stp, verbose = verbose, kwargs...)
end

function TR_Nwt_L(h, a, b, α, β, stp; verbose = false, kwargs...)
    return TR_N_L(h, a, stp, verbose = verbose, kwargs...)
end

function TR_cub2_L(h, a, b, α, β, stp; verbose = false, kwargs...)
    return TR_cub2_L(h, a, stp, verbose = verbose, kwargs...)
end

function TR_seccub2_L(h, a, b, α, β, stp; verbose = false, kwargs...)
    return TR_seccub2_L(h, a, stp, verbose = verbose, kwargs...)
end

