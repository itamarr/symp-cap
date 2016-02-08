function [Fun, gradFun] = FuncToMinimize(coeffs,P,m,n)

    Fun = F(coeffs,P,m,n);
    gradFun = dF(coeffs,P,m,n);
end
