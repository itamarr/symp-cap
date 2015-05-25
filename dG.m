function [y] = dG(x,n)
%DG Describes an ellipsoid

a=[1;2;3;4;5;6;7;8;9;10];
if(size(a,1)~=n)
    error('function G: Dimension mismatch!');
end

y=zeros(2*n,1);
for i=1:n
    y(i)=(a(i)^2)*x(i);
    y(n+i)=(a(i)^2)*x(n+i);
end

y=y/2;
 %Lp-norm
 %p=3;
 %y=((abs(x)/norm(x,p)).^(p-2)).*x/2;

end
