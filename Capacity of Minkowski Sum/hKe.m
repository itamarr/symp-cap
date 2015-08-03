function [ y ] = hKe( x, P, n, e)
%GMinkSum Computes the h_{K_e} function for the (K polar + B(e))-polar body.
%       The formula is (||x||_std * ||x||_{K polar} ) / (||x||_std +
%       e*||x||_{K polar})
%
%       Parameters:
%       x - The point whose norm we want to compute.
%       P - The polytope given like so: Each column is a dimension, each
%       row is an extreme point of the polytope.
%       n - Half the dimension
%       e - The radius of the ball B(e) in the formula above.

nx = norm(x);
normKp = NormOfPolarBody(P,x);
nom = nx * normKp;
denom = norm(x) + e*nx;

y = nom / denom;


end

