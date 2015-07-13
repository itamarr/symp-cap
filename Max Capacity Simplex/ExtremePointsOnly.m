function [ extP ] = ExtremePointsOnly( P )
%EXTREMEPOINTSONLY Gets a list of points as row vectors of a matrix,
%returns a list of the extreme points of the convex hull.
convP = convhulln(P);
X = unique(reshape(convP,1,[])); %get the list of unique row indices appearing in the faces of convex hull.

extP = P(X,:); %select rows according to index;

end

