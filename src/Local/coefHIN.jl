function coefHIN3(x,dy)         # Cubic Hermite interpolation (Newton)
    # vector x: [x_1,x_2]
    # matrix dy contains 2 values, function and derivative each of the 2 points 
    n, m = size(dy)              
    dd = DividedDifference2(x,dy)   # get the divided difference array

    α = x[2] - x[1]
    a = zeros(typeof(x[1]),4)
    a[1] = dd[1,1]
    a[2] = dd[2,2]
    a[3] = dd[3,3] - dd[4,4]*α
    a[4] = dd[4,4]
    
    return a
end


function coefHIN5(x,dy)         # Quintic Hermite interpolation (Newton)
    # vector x: [x_1,x_2]
    # matrix dy contains 3 values, function, derivative and second derivatives of each of the 2 points 
    n, m = size(dy)              
    dd = DividedDifference2(x,dy)   # get the divided difference array

    α = x[2] - x[1]
    a = zeros(typeof(x[1]),6)
    a[1] = dd[1,1]
    a[2] = dd[2,2]
    a[3] = dd[3,3]
    a[4] = dd[4,4] - α*dd[5,5] + α^2*dd[6,6]
    a[5] = dd[5,5] - 2*α*dd[6,6]
    a[6] = dd[6,6]
    
    return a
end


