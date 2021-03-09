function DividedDifference2(x,dy)  # generate array of divided differences
    n, m = size(dy)                # n data points, m derivatives (0 to m-1)
    dd = zeros(typeof(x[1]),n*m,n*m)               # matrix of divided differences
    z = zeros(typeof(x[1]),n*m) 
    k=1;
    for i = 1:n                    # n data points
        for j = 1:m                # m derivatives (0 to m-1) at each point
            k = (i-1)*m + j        # row index
            z[k] = x[i]      
            dd[k,1] = dy[i,1]      # 0th divided difference in first column
            for l = 2:k            # column index for the remaining columns
                if dd[k,l-1] == dd[k-1,l-1]  # left and top-left neighbors are repeated
                    if l<=m
                        dd[k,l] = dy[i,l] / factorial(l-1)
                    else dd[k,l] = big(0.0)
                    end
                else                       
                    dd[k,l] = (dd[k,l-1] - dd[k-1,l-1]) / (z[k]-z[k-l+1])
                end
            end
        end
    end
    return dd
end

