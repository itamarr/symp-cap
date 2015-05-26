function [y] = dG(x,n,p)
%DG Summary of this function goes here
%   Detailed explanation goes here
%y=[x(1)^3;x(2)^3]/(2*sqrt(x(1)^4+x(2)^4));
y=((abs(x)/norm(x,p)).^(p-2)).*x/2;
end
