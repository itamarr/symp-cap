function [y] = dF(x,P,m,n)
%DF This function is the gradient of the action function, denoted by F in
%Going-Jaeschke's paper.
%   Parameters:
%   x - A vector with size [2*n*m 1], such that each m entries represent a
%   vector in the space of derivatives of curves. This represents the
%   (piecewise) derivatives of a closed characteristic of the body P.
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix where k is the number of vertices
%   in the convex hull.
%   m - The number of subintervals of [0,1], as in the paper.
%   n - The half dimension of the body represented in P.

J2n = [zeros(n),eye(n);-eye(n),zeros(n)];

y = zeros(2*n*m,1);
for j = 0:m-1
    y(j*2*n+1:(j+1)*2*n)=J2n*dG(-J2n*x(j*2*n+1:(j+1)*2*n),P,n);
end

y = y/m;

end

