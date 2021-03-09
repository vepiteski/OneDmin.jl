#abstract type TR1DModel{T} end

abstract type CubicModel{T} <: TR1DModel{T} end

abstract type CubicModel2Pt{T} <: CubicModel{T} end

"Cubic model solution"
function solveModel(m:: CubicModel{T}, t :: T) where T
    disc = m.H^2 - T(2)*m.C*m.g
    d = disc < 0 ? -sign(m.C)*Inf : (-m.H + sqrt(disc)) / m.C
    return d
end

"Cubic model evaluation"
function m(M:: CubicModel{T}, t :: T) where T
    # Horner rule
    return M.f + t*(M.g + t*( T(0.5)*M.H + t*( T(1.0/6.0)*M.C)))
end

"Cubic model derivative evaluation"
function dm(M:: CubicModel{T}, t :: T) where T
    return M.g  +  t*(M.H  +  t*T(0.5)*M.C)
end

"Cubic model predicted reduction"
function Δm(M:: CubicModel{T}, t :: T) where T
    return M.g * t  +  T(0.5) * M.H * t^2  + T(1.0/6.0)*M.C*t^3
end

"Cubic model predicted reduction given two points"
function Δm(M:: CubicModel{T}, t1 :: T, t2 :: T) where T
    c1 = M.g * t1  +  T(0.5) * M.H * t1^2  + T(1.0/6.0)*M.C*t1^3
    c2 = M.g * t2  +  T(0.5) * M.H * t2^2  + T(1.0/6.0)*M.C*t2^3
    return c1 - c2
end

"Taylor cubic polynomial model"
mutable struct Taylor3{T} <: CubicModel{T}
    f :: T
    g :: T
    H :: T
    C :: T
end
Taylor3{T}() where T = Taylor3{T}(0,0,0,0)


"2 points interpolation cubic polynomial model"
mutable struct TwoPoints3{T} <: CubicModel2Pt{T}
    f :: T
    g :: T
    H :: T
    C :: T
end
TwoPoints3{T}() where T = TwoPoints3{T}(0,0,0,0)



function initModel(M:: Taylor3{T}, ϕ :: LineModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g = dϕt
    M.H = hess(ϕ, t)
    M.C = tens(ϕ,t)
    return M
end

function updateModel!(M:: Taylor3{T}, ϕ :: LineModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g  = dϕt
    M.H = hess(ϕ, t)
    M.C = tens(ϕ,t)
    return M
end


"Secant cubic model"
mutable struct secant3{T} <: CubicModel2Pt{T}
    f :: T
    g :: T
    H :: T
    C :: T
    Q :: T  # to test plain secant
end
secant3{T}() where T = secant3{T}(0,0,0,0,0)


function initModel(M:: CubicModel2Pt{T}, ϕ :: LineModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g = dϕt
    M.H, M.C = initHC_s!(ϕ, t, ϕt, dϕt)
    return M
end

function updateModel!(M:: CubicModel2Pt{T}, ϕ :: LineModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g = dϕt
    M.H, M.C = updateHC_s!(ϕ, t, ϕt, dϕt)
    return M
end

function initModel(M:: secant3{T}, ϕ :: LineModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g = dϕt
    M.H, M.C, M.Q = initHC_s!(ϕ, t, ϕt, dϕt)
    return M
end

function updateModel!(M:: secant3{T}, ϕ :: LineModel, t :: T, ϕt :: T, dϕt :: T) where T
    M.f = ϕt
    M.g = dϕt
    M.H, M.C, M.Q = updateHC_s!(ϕ, t, ϕt, dϕt)
    return M
end


let (tpred, ϕpred, dϕpred) = (NaN, NaN, NaN)
    
    global function initHC_s!(ϕ :: LineModel, t :: T, ϕt :: T, dϕt :: T) where T
        (tpred, ϕpred, dϕpred) = (t, ϕt, dϕt)
        return T(1), T(0), T(1)
    end
       
    global function updateHC_s!(ϕ :: LineModel, t :: T, ϕt :: T, dϕt :: T) where T
        # estimate quadratic and cubic coefficients from the hermite interpolant
        a = coefHIN3([t;tpred], [ϕt dϕt; ϕpred dϕpred])
        s = t - tpred
        y = dϕt - dϕpred
        Q = y/s

        tpred, ϕpred, dϕpred = t, ϕt, dϕt
        
        return T(2.0)*a[3], T(6.0)*a[4], Q
    end

    "Two consecutive secant steps on the cubic model"
    global function solveModel(M:: secant3{T}, t :: T) where T
        # the "Taylor" cubic polynomial fits with the function and its derivative at t and tpred
        # Thus, secant iteration may be performed on the cubic model without knowledge of ϕ.
        # first secant iteration
        #s = t - tpred
        dmt = M.g  #  dm(M, t)
        #dmpred = dm(M, tpred)
        #y = dmt - dmpred
        dS = - (dmt / M.Q)

        # second secant iteration
        t2 = t + dS
        dm2 = dm(M, t2)
        y2 = dm2 - dmt
        d2 = - (dm2 * dS / y2)
        dS2 = dS + d2
        
        #return dS # for testing, should be equivalent to secant
        #return dS2 # for testing, should be equivalent to secant

        # try the cubic's root
        disc = M.H^2 - 2.0*M.C*M.g
        dC = disc < 0 ? -sign(M.C)*Inf : (-M.H + sqrt(disc)) / M.C
        return dC
    end
end
