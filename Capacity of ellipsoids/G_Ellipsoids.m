function [ y ] = G_Ellipsoids(x,R,n)
%G_Ellipsoids Calculates the norm of the point represented by x with
%respect to the ellipsoid given by R.

if (size(R,2) ~= n)
    error('function G_Ellipsoids: Dimension mismatch!');
end

inverseRadiiSquared = R.^(-2);

% This calculation assumes the coordinates are q_1,p_1,...,q_n,p_n
%r = repmat(inverseRadiiSquared, n, 1);
%y = 0.25 * r(:)' * x.^2;

% This calculation assumes the coordinates are q_1,..,q_n,p_1,..,p_n
r = repmat(inverseRadiiSquared, 1, 2);
y = 0.25 * r * x.^2;


%sum = 0;
%for i = 1:2:2*n
%    sum = sum + inverseRadiiSquared(i) * (x(i)^2 + x(i+1)^2);
%end

%y = 0.25 * sum;
