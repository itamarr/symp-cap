function [y] = F(x,m,n)
%F Summary of this function goes here
%   Detailed explanation goes here

mJ2n = [zeros(n),eye(n);-eye(n),zeros(n)];

sum = 0;
for j=0:m-1
    sum=sum+G(mJ2n*x(j*2*n+1:(j+1)*2*n),n);
end

y=sum/m^2;

end
