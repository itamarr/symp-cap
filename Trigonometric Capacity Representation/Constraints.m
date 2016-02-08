function [c, ceq, gradc, gradceq] = Constraints(coeffs,m,n)
%Constraints This function returns the constrains required for fmincon used
%in Capacity function.
    global constraintMatrix;
    
    c = [];
    gradc = [];
    
    ceq = coeffs' * constraintMatrix * coeffs - 1;
    gradceq = (constraintMatrix  + constraintMatrix')*coeffs;
    
end
