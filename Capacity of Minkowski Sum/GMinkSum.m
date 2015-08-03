function [ y ] = GMinkSum( x, P, n, e )
%GMinkSum Calculates the G function for hKe.
%       For posterity: MinkSum is supposed to relate to Minkowsky sum of
%       the polar body with a ball. The goal is to eliminate flat faces in
%       the original body.

h = hKe(x, P, n, e);
y = 0.25 * h^2;

end

