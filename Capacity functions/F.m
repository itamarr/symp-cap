function [y] = F(x,P,m,n)
%F Calculates the value of the functional F as in Going-Jaweschke's paper.
%This is basically the action of the characteristic whose derivative is
%represented in x, up to a scalar multiple.
%   Parameters:
%   x - A vector with size [2*n*m 1], such that each m entries represent a
%   vector in the space of derivatives of curves. This represents the
%   (piecewise) derivatives of a closed characteristic of the body P.
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix where k is the number of vertices
%   in the convex hull.
%   m - The number of subintervals of [0,1], as in the paper.
%   n - The half dimension of the body represented in P.

%global mJ2n
%if (isempty(mJ2n))
	mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)];
%end

sum = 0;
for j=0:m-1
    sum=sum+G(mJ2n*x(j*2*n+1:(j+1)*2*n),P,n);
end
y=sum/m;



end

