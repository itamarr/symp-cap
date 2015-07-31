function [y] = dF(x,P,m,n)
%DF Summary of this function goes here
%   Detailed explanation goes here

%global J2n
%if (isempty(J2n) || size(J2n,1) ~= 2*n)
    J2n = [zeros(n),eye(n);-eye(n),zeros(n)];
%end

y=zeros(2*n*m,1);
for j=0:m-1
    y(j*2*n+1:(j+1)*2*n)=J2n*dG(-J2n*x(j*2*n+1:(j+1)*2*n),P,n);
end

y=y/m;

end

