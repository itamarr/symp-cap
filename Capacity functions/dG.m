function [y] = dG(x,P,n)
%DG The gradient of the square of the norm of the polar body to the body represented as
%P. The function of which this is a gradient is denoted by G in the paper
%by Going-Jaeschke.
%   Parameters:
%   x - A vector with size [2*n 1] in the dual space to where the body P
%   lies.
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix where k is the number of vertices
%   in the convex hull.
%   n - The half dimension of the body represented in P.

y = (1/2) * NormOfPolarBody(P, x) * (GradOfPolarNorm(P, x)');

end
