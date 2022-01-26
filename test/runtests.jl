using OneDmin
using Test
using Logging

fh(x) = - ( x + sin(x)) * exp(-x^2)
x0 = big(0.2); x1 = big(1.0)

using ForwardDiff

gh(x) = ForwardDiff.derivative(fh, x)  # first derivative
Hh(x) = ForwardDiff.derivative(gh, x)  # second derivative
Ch(x) = ForwardDiff.derivative(Hh, x)  # third derivative

setprecision(4096)

function test(algo::Function, x, xm)
    h = fh(x)
    dh = gh(x)
    d2h = Hh(x)
    d3h = Ch(x)
    
    hm = fh(xm)
    dhm = gh(xm)
    d2hm = Hh(xm)
    
    it = 0
    maxit = 25
    #with_logger(Logging.NullLogger()) do
    #    @info "it  x  |h'|"
        while abs(dh) > BigFloat("1.0e-1200")
            xp = algo(x,xm,h,hm,dh,dhm,d2h,d2hm,d3h,0)
            x, xm = xp, x
            
            hm = h
            dhm = dh
            d2hm = d2h
            
            h = fh(x)
            dh = gh(x)
            d2h = Hh(x)
            d3h = Ch(x)
            
            it += 1
            if (it > maxit)  break  end
    #        @info @sprintf("%2d",it), @sprintf("%f",x), @sprintf("%.2e",abs(dh))
        end
    #end
    return it 
    
end


iter1 = Dict(base .=> [16,12,8])
iter2 = Dict(quintic .=> [8,7])
iter3 = Dict(cubic1 .=> [8,8,11])
iter4 = Dict(cubic0 .=> [11,11,15])

iter = merge(iter1, iter2, iter3, iter4)

@testset "OneDmin.jl" begin
    xm = x1; x = x0
    
    @info "Testing algos in base."
    for algo in base
        @info "testing" algo
        it = test(algo, x, xm)
        @test it == iter[algo]
    end
    
    @info "Testing algos in quintic."
    for algo in quintic
        @info "testing" algo
        it = test(algo, x, xm)
        @test it == iter[algo]
    end
    
    @info "Testing algos in cubic1."
    for algo in cubic1
        @info "testing" algo
        it = test(algo, x, xm)
        @test it == iter[algo]
    end
    
    @info "Testing algos in cubic0."
    for algo in cubic0
        @info "testing" algo
        it = test(algo, x, xm)
        @test it == iter[algo]
    end

    include("test1D.jl")
    
    
end
