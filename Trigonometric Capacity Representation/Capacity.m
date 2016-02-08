function [c, gamma] = Capacity(P,n,m)

    C = centroid(P);
    P = P - repmat(C,size(P,1),1);
    
    global constraintMatrix;
    sz = 2*n*(2*m + 1);
    constraintMatrix = zeros(sz,sz);
    
    for k = 1:m
        for j = 1:n
            if (((n+j-1)*(2*m+1)+1+k > sz) || ((j-1)*(2*m+1)+1+m+k > sz) || ((j-1)*(2*m+1)+1+k > sz) || ((n+j-1)*(2*m+1)+1+m+k > sz))
                print "wtf"
            end
            constraintMatrix((n+j-1)*(2*m+1)+1+k, (j-1)*(2*m+1)+1+m+k) = k;
            constraintMatrix((j-1)*(2*m+1)+1+k, (n+j-1)*(2*m+1)+1+m+k) = -k;
        end
    end
    
    % generating starting path's coefficients.
    
    coeffs = (1/100)*rand(2*n*(2*m+1),1);
    coeffs(1)=0;
    coeffs(2)=1;
    coeffs((2*m+1)+1+m+1)=-1;
    
    %coeffs = rand(2*n*(2*m+1),1);
    l = coeffs'*constraintMatrix*coeffs;
    while l<=0
        print('beep');
        coeffs = rand(2*n*(2*m+1),1);
        l = abs(coeffs'*constraintMatrix*coeffs);
    end
    coeffs = coeffs / sqrt(l);

    
    options = optimoptions('fmincon','GradObj','on','GradConstr','on');
    options.Display = 'iter';
    options.Algorithm = 'interior-point'; %% should try which works best, Maybe this is better that sqp?
    options.MaxIter = 300;
    options.TolCon = 1e-7;
    options.TolFun = 1e-7;
    options.TolX = 1e-7;
    
    cond = @(x) Constraints(x,m,n);
    pf = @(x) FuncToMinimize(x,P,m,n);

    [g, udotval] = fmincon(pf,coeffs,[],[],[],[],[],[],cond,options);
    
    % this is the action of the minimizer
    c = 2*F(g,P,m,n);
    gamma = g;
end