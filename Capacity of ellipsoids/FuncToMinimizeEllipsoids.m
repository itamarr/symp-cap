function [Fun, gradFun] = FuncToMinimizeEllipsoids(x,R,m,n)
%
    %avg = P'*ones(2*n+1,1) / (2*n+1);
    %P = P - repmat(avg',2*n+1,1);
    Fun = F_Ellipsoids(x,R,m,n);
    gradFun = dF_Ellipsoids(x,R,m,n);
end