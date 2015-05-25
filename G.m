function [ y ] = G(x,n)
% %G Describes an ellipsoid
% %   
% a=[1;2;3;4;5;6;7;8;9;10];
% if(size(a,1)~=n)
%     error('function G: Dimension mismatch!');
% end
% 
% sum=0;
% for i=1:n
%     sum=sum+(a(i)^2)*(x(i)^2+x(i+n)^2);
% end
% 
% y=sum/4;
% 
% %Lp-norm
% %p=3;
%y=(norm(x,2)^2)/4;
% 
%%l1ball = [1 0;-1 0;0 1;0 -1];

%cube = [2 2; -2 2; -2 -2; 2 -2];
std2splx = [1/2 -1/2; -1/2 1/2; -1/2 -1/2];
y = (NormOfPolarBody(std2splx, x')^2)/4;

end
