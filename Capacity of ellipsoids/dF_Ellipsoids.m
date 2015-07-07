function [y] = dF_Ellipsoids(x,R,m,n)
%DF Summary of this function goes here
%   Detailed explanation goes here

J2n = [zeros(n),eye(n);-eye(n),zeros(n)];

y=zeros(2*n*m,1);
for j=0:m-1
    y(j*2*n+1:(j+1)*2*n) = J2n*dG_Ellipsoids(-J2n*x(j*2*n+1:(j+1)*2*n),R,n);
end

y = y/m;

end

