function [ grad ] = GradOfPolarNorm( P, y )
%GradOfPolarNorm Calculates the gradient of the norm of polar body to P at
%point y.
%   Parameters:
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix where k is the number of vertices
%   in the convex hull.
%   y - The point at which we want to calculate the gradient. This is
%   represented as a column vector.

    [M,I] = max(P*y); %vertices of P correspond to normal vectors of P^o,
                      %Find which has smallest angle w.r.t y.
    grad = P(I(1),:);
  
end
