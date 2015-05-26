function [ y ] = G(x,n,p)
%G Summary of this function goes here
%  
%y=sqrt(x(1)^4+x(2)^4)/4;
y=norm(x,p)^2/4;
end
