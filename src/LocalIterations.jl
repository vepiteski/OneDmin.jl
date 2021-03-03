# A bunch of functions to perform a single iteration to minimize h(t)
# and/or compute a root of h'(t)
include("DividedDifference2.jl")
include("coefHIN.jl")


# basic iterations
base =  Function[]


"""One iteration of the secant method"""
function s_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    s = xm - x
    y = dhm - dh

    d = - (dh * s / y)

    xp = x + d
    return xp
end
push!(base, s_iteration)


"""One iteration of the Newton's method"""
function N_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    d = - (dh / d2h)

    xp = x + d
    return xp
end
push!(base, N_iteration)


"""One iteration of a cubic Taylor approximation
   may also be described as a quadratic Taylor approximation of h'
"""
function C_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,d3h,d3hm)
    # solution cubique...
    disc = d2h^2 - 2.0*d3h*dh
    d = disc < 0 ? -sign(d3h)*Inf : (-d2h + sqrt(disc)) / d3h

    xp = x + d
    return xp
end
push!(base, C_iteration)



# iterations on the quintic interpolant of h, h' and h''
quintic = Function[]





"""Two Newton iteration on the quintic interpolant of h
   at the current and previous point
"""
function N2Opt_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    d = - (dh / d2h)
    x2 = x + d
    # estimate dh2 from the hermite interpolant

    a = coefHIN5([x;xm],[h dh d2h; hm dhm d2hm])

    p = Polynomial(a)
    dp = Polynomials.derivative(p)
    d2p = Polynomials.derivative(dp)

    d2 = -(dp(d) / d2p(d))

    xp = x2 + d2

    return xp
end
push!(quintic, N2Opt_iteration)


"""Three Newton iteration on the quintic interpolant of h
   at the current and previous point
"""
function N3Opt_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    d = - (dh / d2h)
    x2 = x + d
    # estimate dh2 from the hermite interpolant

    a = coefHIN5([x;xm],[h dh d2h; hm dhm d2hm])

    p = Polynomial(a)
    dp = Polynomials.derivative(p)
    d2p = Polynomials.derivative(dp)

    d2 = -(dp(d) / d2p(d))

    x3 = x2 + d2
    d3 = -(dp(d+d2) / d2p(d+d2))

    xp = x3 + d3
    return xp
end
push!(quintic, N3Opt_iteration)


# iterations on the cubic interpolant of h' and h''
cubic1 = Function[]


"""Two Newton iteration on the cubic interpolant of h'
   at the current and previous point
"""
function N2_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)

    d = - (dh / d2h)
    x2 = x + d
    # estimate dh2 from the hermite interpolant

    a = coefHIN3([x;xm],[dh d2h; dhm d2hm])
    p = Polynomial(a)
    dp = Polynomials.derivative(p)

    d2 = -(p(d) / dp(d))

    xp = x2 + d2
    return xp
end
push!(cubic1, N2_iteration)


"""Three Newton iterations on the cubic interpolant of h'
   at the current and previous point
"""
function N3_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    d = - (dh / d2h)
    x2 = x + d
    # estimate dh2 from the hermite interpolant

    a = coefHIN3([x;xm],[dh d2h; dhm d2hm])

    p = Polynomial(a)
    dp = Polynomials.derivative(p)

    d2 = -(p(d) / dp(d))

    x3 = x2 + d2
    d3 = -(p(d+d2) / dp(d+d2))
    xp = x3 + d3
        
    return xp
end
push!(cubic1, N3_iteration)


"""Two secant iterations on the cubic interpolant of h'
   at the current and previous point
"""
function s2Opt_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    s = x - xm
    y = dh - dhm

    d = - (dh * s / y)

    x2 = x + d
    # estimate dh2 from the hermite interpolant
    a = coefHIN3([x;xm],[h dh; hm dhm])
    s2 = x2 - x

    p = Polynomial(a)
    dp = Polynomials.derivative(p)
    dh2 = dp( s2)

    # perform a secant iteration with this estimated dh2
    y2 = dh2 - dh
    d2 = - (dh2 * s2 / y2)

    xp = x2 + d2
    if ( xp == x ) xp = x2 end #safeguard
    return xp
end
push!(cubic1, s2Opt_iteration)



# iterations on the cubic interpolant of h and h'
cubic0 = Function[]




"""Minimizer of the cubic interpolant of h and dh"""
function sCOpt_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    # Hermite cubic interpolant
    a = coefHIN3([x;xm],[h dh; hm dhm])
    p = Polynomial(a)
    dp = Polynomials.derivative(p)
    dh  = a[2]
    d2h = 2.0*a[3]
    d3h = 3.0*a[4]
    disc = d2h^2 - 2.0*d3h*dh
    d = disc < 0 ? -sign(d3h)*Inf : (-d2h + sqrt(disc)) / d3h

    xp = x + d

    return xp
end
push!(cubic0, sCOpt_iteration)


"""Two Newton's iterations on the cubic interpolant of h
   at the current and previous point
"""
function sNOpt_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    # Hermite cubic interpolant
    a = coefHIN3([x;xm],[h dh; hm dhm])
    p = Polynomial(a)
    dp = Polynomials.derivative(p)
    d2p = Polynomials.derivative(dp)

    # Newton step on the interpolant
    # Do not forget that the interpolant is rooted at 0
    d = - dp(0.0) / d2p(0.0)
    d2 = - dp(d) / d2p(d)
    xp = x + d + d2
    
    return xp
end
push!(cubic0, sNOpt_iteration)


"""Three secant iterations on the cubic interpolant of h
   at the current and previous point; not efficient????
"""
function s3Opt_iteration(x,xm,h,hm,dh,dhm,d2h,d2hm,::Any,::Any)
    s = x - xm
    y = dh - dhm

    d = - (dh * s / y)

    x2 = x + d
    
    # estimate dh2 from the hermite interpolant
    a = coefHIN3([x;xm],[h dh; hm dhm])
    p = Polynomial(a)
    dp = Polynomials.derivative(p)

    s2 = x2 - x
    dh2 = dp(s2)
    # perform a secant iteration with this estimated dh2
    y2 = dh2 - dh
    d2 = - (dh2 * s2 / y2)
    x3 = x2 + d2

    s3 = x3 - x2
    dh3 = dp(s3)
    # perform a secant iteration with this estimated dh3
    y3 = dh3 - dh2
    d3 = - (dh3 * s3 / y3)

    xp = x3 + d3
    # not good, seems to fall back close to regular secant iteration
    if ( xp == x ) xp = x3 end#(x+xm)/2.0 end#@show a end
    return xp
end
push!(cubic0, s3Opt_iteration)


export base
for n in nameof.(base)
    @eval export $n
end

export quintic
for n in nameof.(quintic)
    @eval export $n
end

export cubic1
for n in nameof.(cubic1)
    @eval export $n
end

export cubic0
for n in nameof.(cubic0)
    @eval export $n
end

