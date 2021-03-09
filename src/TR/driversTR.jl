

function TR1D_N(h, a, b, α, β, stp; kwargs...)
    return TR1D(h, a=0.0, b=Inf, α=α, β=β, stp=stp; kwargs...)
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

