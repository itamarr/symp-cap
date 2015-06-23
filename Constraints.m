function [c, ceq, gradc, gradceq] = Constraints(x,m,n,A2n)
    c=[];
    gradc=[];
    
    ceq = ((x'*A2n*x)/m^2) - 1;
    gradceq = (1/m^2)*(A2n+A2n')*x;
    
end