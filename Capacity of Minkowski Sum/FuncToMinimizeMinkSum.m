function [Fun, gradFun] = FuncToMinimizeMinkSum(x,P,m,n,e)
%
    %avg = P'*ones(2*n+1,1) / (2*n+1);
    %P = P - repmat(avg',2*n+1,1);
    Fun = FMinkSum(x,P,m,n,e);
    gradFun = dFMinkSum(x,P,m,n,e);
end
