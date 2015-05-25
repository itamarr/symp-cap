function [ y ] = G(x,n)
%G Describes an ellipsoid
%   
a=[1;2;3;4;5;6;7;8;9;10];
if(size(a,1)~=n)
    error('function G: Dimension mismatch!');
end

sum=0;
for i=1:n
    sum=sum+(a(i)^2)*(x(i)^2+x(i+n)^2);
end

y=sum/4;

%Lp-norm
%p=3;
%y=(norm(x,p)^2)/4;


end
