function [ y ] = dGMinkSum( x, P, n, e )
%dGMinkSum Computes the gradient of GMinkSum.
%
%       Parameters:
%       x - The point whose norm we want to compute.
%       P - The polytope given like so: Each column is a dimension, each
%       row is an extreme point of the polytope.
%       n - Half the dimension
%       e - The radius of the ball B(e) in the formula above.

y = 0.5 * hKe(x, P, n, e) * dhKe(x, P, n, e);

end

