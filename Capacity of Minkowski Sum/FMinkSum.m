function [ y ] = FMinkSum( x,P,m,n,e )
%FMinkSum Summary of this function goes here
%   Detailed explanation goes here

%global mJ2n
%if (isempty(mJ2n))
	mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)];
%end

sum = 0;
for j=0:m-1
    sum=sum+GMinkSum(mJ2n*x(j*2*n+1:(j+1)*2*n),P,n,e);
end
y=sum/m;



end

