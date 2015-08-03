function [y] = dG(x,P,n)
%DG Describes an ellipsoid
% 

y = (1/2)*NormOfPolarBody(P, x)*(GradOfPolarNorm(P,x)');

end
