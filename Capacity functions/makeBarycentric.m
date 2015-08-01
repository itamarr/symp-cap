function [barycentricP] = makeBarycentric(P)
%MAKEBARYCENTRIC makes a polytope barycentric
%   vertices should be written in rows
s=size(P,1);
barycenter=sum(P,1)/s;
barycentricP=P-repmat(barycenter,s,1);

end
