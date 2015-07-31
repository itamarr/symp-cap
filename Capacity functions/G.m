function [ y ] = G(x,P,n)
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
%l1ball = [1 0;-1 0;0 1;0 -1];

% rad=1;
% N=4;
% theta = [0 : 2*pi/N : 2*pi];
% Pos = rad * exp(1i*theta);
% X = real(Pos);
% Y = imag(Pos);
% regpoly = [X ; Y]';
%cube = [1 1; -1 1; -1 -1; 1 -1];
%std2splx = [1/2 -1/2; -1/2 1/2; -1/2 -1/2];
%stdSimplex = [zeros(1,4);eye(4)] - 1/8;  %% 1/4 1/4 ... 1/4 is on the face, take half distance along this line

y = (NormOfPolarBody(P, x)^2)/4;

end
