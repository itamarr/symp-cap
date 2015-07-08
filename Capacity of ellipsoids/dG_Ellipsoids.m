function [y] = dG_Ellipsoids(x,R,n)
%dG_Ellipsoids Calculates the gradient of the norm of the point represented by x with
%respect to the ellipsoid given by R.

if (size(R,2) ~= n)
    error('function G_Ellipsoids: Dimension mismatch!');
end

radiiSquared = R.^2;

% This calculation assumes the coordinates are q_1,p_1,...,q_n,p_n
%r = repmat(inverseRadiiSquared, n, 1);
%y = 0.5 * (r(:).* x.^2);

% This calculation assumes the coordinates are q_1,..,q_n,p_1,..,p_n
r = repmat(radiiSquared, 1, 2);
y = 0.5 * (r'.* x); % elementwise r_i times x_i


%y = zeros(2*n,1);
%for i = 1:2:2*n
%    y(i) = inverseRadiiSqaured(i) * x(i);
%    y(i+1) = inverseRadiiSqaured(i) * x(i+1);
%end

%y = 0.5 * y;
