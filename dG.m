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
% y=x/2;
%  %Lp-norm
%p=2;
%y=((abs(x)/norm(x,p)).^(p-2)).*x/2;
%l1ball = [1 0;-1 0;0 1;0 -1];
rad=1;
N=4;
theta = [0 : 2*pi/N : 2*pi];
Pos = rad * exp(1i*theta);
X = real(Pos);
Y = imag(Pos);
regpoly = [X ; Y]';
%cube = [1 1; -1 1; -1 -1; 1 -1];
%std2splx = [1/2 -1/2; -1/2 1/2; -1/2 -1/2];
stdSimplex = [zeros(1,4);eye(4)] - 1/8;
y = (1/2)*NormOfPolarBody(stdSimplex, x')*(GradOfPolarNorm(stdSimplex,x)');

end
