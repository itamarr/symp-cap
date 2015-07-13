function [ centP ] = CentralizePolytope( P )
%CENTRALIZESIMPLEX Gets a simplex as list of row vectors of extreme points,
%translates it such that the center of mass is at zero.
%   Detailed explanation goes here
C = centroid(P);
centP = P - repmat(C,size(P,1),1);
end

