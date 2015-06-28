function [Fun, gradFun] = FuncToMinimize(x,P,m,n)
%
    %avg = P'*ones(2*n+1,1) / (2*n+1);
    %P = P - repmat(avg',2*n+1,1);
    Fun = F(x,P,m,n);
    gradFun = dF(x,P,m,n);
end