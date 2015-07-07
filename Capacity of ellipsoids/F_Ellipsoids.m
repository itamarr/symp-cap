function [y] = F_Ellipsoids(x,R,m,n)
%F Summary of this function goes here
%   Detailed explanation goes here

mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)]; %'mJ2n' is -J_2n.

sum = 0;
for j=0:m-1
    sum = sum + G_Ellipsoids(mJ2n*x(j*2*n+1:(j+1)*2*n),R,n);
end
y = sum/m;



end

