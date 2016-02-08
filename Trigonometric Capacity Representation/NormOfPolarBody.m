function [ n ] = NormOfPolarBody( P, y )
%NormOfPolarBody Calculates the norm of y in relation to the polar of P.
%   Parameters:
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix where k is the number of vertices
%   in the convex hull.
%   y -  The point at which we want to calculate the gradient. This is
%   represented as a column vector.


    [M,I] = max(P*y);
    
    n = (P(I(1),:) * y);
end

