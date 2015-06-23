function [Fun, gradFun] = FuncToMinimize(x,m,n)
    Fun = F(x,m,n);
    gradFun = dF(x,m,n);
end