function [ y ] = G(x,P,n)
% %G Describes an ellipsoid
% %   

y = (NormOfPolarBody(P, x)^2)/4;

end
