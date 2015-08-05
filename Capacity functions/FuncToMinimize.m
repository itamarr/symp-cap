function [Fun, gradFun] = FuncToMinimize(x,P,m,n)
% FuncToMinimize This function plugs the arguments into the F and dF
% functions. This is then passed as a parameter to fmincon.
%   Parameters:
%   x - A vector with size [2*n*m 1], such that each m entries represent a
%   vector in the space of derivatives of curves. This represents the
%   (piecewise) derivatives of a closed characteristic of the body P.
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix where k is the number of vertices
%   in the convex hull.
%   m - The number of subintervals of [0,1], as in the paper.
%   n - The half dimension of the body represented in P.
 
    Fun = F(x,P,m,n);
    gradFun = dF(x,P,m,n);
end
