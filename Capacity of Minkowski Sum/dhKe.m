function [ y ] = dhKe( x, P, n, e )
%dhKe The grad function of hKe.
%       The formula is (e||x||_{K polar} * \grad ||x||_std + ||x||_std^2*(\grad ||x||_{K polar}) /
%       (||x||_std + e*||x||_{K polar})^2

nK = NormOfPolarBody(P, x);
nx = norm(x);
dNormx = x / nx;
dNormKx = GradOfPolarNorm(P, x)';

y = (e*nK^2*dNormx + nx^2*dNormKx) / (nx + e*nK)^2;



end

