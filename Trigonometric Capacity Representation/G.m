function [ y ] = G(x,P,n)
%G The square of the norm of the polar body to P. This is denoted by G as
%in the paper by Going-Jaeschke.
%   Parameters:
%   x - A vector with size [2*n 1] in the dual space to where the body P
%   lies.
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix where k is the number of vertices
%   in the convex hull.
%   n - The half dimension of the body represented in P.

y= zeros(1, size(x,2));

for i=1:size(x,2)
    y(i) = (NormOfPolarBody(P, x(:,i))^2)/4;
end

end