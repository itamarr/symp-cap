function [y] = dG(x,n)
%DG Describes an ellipsoid
% 
% a=[1;2;3;4;5;6;7;8;9;10];
% if(size(a,1)~=n)
%     error('function G: Dimension mismatch!');
% end
% 
% y=zeros(2*n,1);
% for i=1:n
%     y(i)=(a(i)^2)*x(i);
%     y(n+i)=(a(i)^2)*x(n+i);
% end
% 
% y=y/2;
%  %Lp-norm
%p=2;
%y=((abs(x)/norm(x,p)).^(p-2)).*x/2;
%l1ball = [1 0;-1 0;0 1;0 -1];

%cube = [2 2; -2 2; -2 -2; 2 -2];
std2splx = [1/2 -1/2; -1/2 1/2; -1/2 -1/2];
y = (1/2)*NormOfPolarBody(std2splx, x')*(GradOfPolarNorm(std2splx,x)');

end
