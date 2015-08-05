function [ centP ] = CentralizePolytope( P )
%CENTRALIZESIMPLEX Translates a polytope such that the center of mass moves
%to zero. 
%   Parameters:
%   P - The extreme points of the poyltope given in a form of a matrix,
%   where the rows of the matrix are the vertices of the polytope.
C = centroid(P);
centP = P - repmat(C,size(P,1),1);
end

