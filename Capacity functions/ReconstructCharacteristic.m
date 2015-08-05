function [ Char ] = ReconstructCharacteristic( x, P, m, n )
%RECONSTRUCTCHARACTERISTIC This function reconstructs the characteristic
%corresponding to the vector of derivatives described in x. This is
%basically done with discrete integration. For reference, see
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

% note that in the paper the elements of (R^2n)^m are indexed as x_0 ...
% x_{m-1}
% A B C stand for each term in formula (2.14)
lambda = F(x,P,m,n);
mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)]; %'mJ2n' is -J_2n.
matX = reshape(x,[2*n,m]);

v = zeros(size(x));


sumOfGrads = 0;
for j = 0:m-1
    sumOfGrads = sumOfGrads + dG(mJ2n*matX(:,j+1),P,n);
end
B = (1/(sqrt(lambda)*m))*sumOfGrads;

tmpV = zeros(2*n,1);
for i = 1:m-1
    for j = 0:i-1
        tmpV = tmpV+matX(:,j+1);
    end
end
C = (-2*sqrt(lambda)/(m^2))*tmpV;

v(1:2*n) = B + C;

for k = 1:m-1
    A = (2*sqrt(lambda)/m)*sum(matX(:,1:k),2); 
    v(k*2*n+1:(k+1)*2*n) = A+B+C;
end

Char = reshape(v,[1,2*n*m]);
end

