using NLPModels
using SolverTools

abstract type TR1DModel{T} end

abstract type QuadraticModel{T} <: TR1DModel{T} end

"Quadratic model solution"
function solveModel(M:: QuadraticModel{T}, t :: T) where T
    return -M.g/M.H
end

"Quadratic model solution in region  Δn ≤ t ≤ Δp "
function solveModel(M:: QuadraticModel{T}, t :: T, Δp :: T, Δn :: T) where T
    dN = -M.g/M.H

    if (Δm(M, Δp, Δn) < 0.0) #| (Δn==0.0)
        d=Δp
        @debug "dp"
    else
        d=Δn
        @debug "dn"
    end
    
    if  (dN < Δp) & (dN > Δn) & (Δm(M, d, dN) > 0.0)
        d=dN 
        @debug "dN = ", d
    end
    
    return d
end


"Quadratic model evaluation"
function m(M:: QuadraticModel{T}, t :: T) where T
    return M.f + M.g * t + T(0.5)*M.H * t^2
end

"Quadratic model predicted reduction"
function Δm(M:: QuadraticModel{T}, t :: T) where T
    return M.g * t  +  T(0.5) * M.H * t^2
end

"Quadratic model predicted reduction given two points"
function Δm(M:: QuadraticModel{T}, t1 :: T, t2 :: T) where T
    #q1 = M.g * t1  +  T(0.5) * M.H * t1^2
    #q2 = M.g * t2  +  T(0.5) * M.H * t2^2
    #diff = q1 - q2
    diff = (M.g + T(0.5) * M.H * (t1+t2)) * (t1 - t2)
    return diff
end

"Taylor quadratic polynomial model"
mutable struct Taylor2{T} <: QuadraticModel{T}
    f :: T
    g :: T
    H :: T
end
Taylor2{T}() where T = Taylor2{T}(0,0,0) 

function initModel(M:: Taylor2{T}, ϕ :: OneDModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g = dϕt
    M.H = hess(ϕ, t)
    return M
end

function updateModel!(M:: Taylor2{T}, ϕ :: OneDModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g  = dϕt
    M.H = hess(ϕ, t)
    return M
end



function initH_N(ϕ :: LineModel, t :: Float64, bidon)
    return hess(ϕ, t)
end

function updateH_N(ϕ :: LineModel, t :: Float64, dϕ :: Float64)
    return hess(ϕ, t)
end




"Secant quadratic model"
mutable struct secant{T} <: QuadraticModel{T}
    f :: T
    g :: T
    H :: T
end
secant{T}() where T = secant{T}(0,0,0) 

function initModel(M:: secant{T}, ϕ :: OneDModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g = dϕt
    M.H = initH_s!(ϕ, t, dϕt)
    return M
end

function updateModel!(M:: secant{T}, ϕ :: OneDModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g = dϕt
    M.H = updateH_s!(ϕ, t, dϕt)
    return M
end


let (tpred, dϕpred) = (NaN, NaN)
    
    global function initH_s!(ϕ :: OneDModel, t :: Float64, dϕ :: Float64)
        (tpred, dϕpred) = (t, dϕ)
        return 1.0
    end
       
    global function updateH_s!(ϕ :: OneDModel, t :: Float64, dϕt :: Float64)
        s = t - tpred
        y = dϕt - dϕpred
        H = y/s

        tpred, dϕpred = t, dϕt
        
        return H
    end
end
