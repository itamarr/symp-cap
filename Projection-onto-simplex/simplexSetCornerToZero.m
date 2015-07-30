function [newSimplexMatrix] = simplexSetCornerToZero(simplexMatrix,indexToZero)
%SIMPLEXSETCORNERTOZERO Translates a simplex to set a corner to zero
%
%   INPUT:
%
%   simplexMatrix = the simplex. Corners in the columns
%
%   indexToZero = corner to be set to zero
%
%   OUTPUT:
%
%   newSimplexMatrix = matrix of the translated simplex, omitting the 

dim = size(simplexMatrix,1);
% if size(simplexMatrix,2) ~= (dim+1)
%     error('Wrong dimensions of simplex!!');
% end

i=indexToZero;

newSimplexMatrix=zeros(dim,size(simplexMatrix,2)-1);
for j=1:i-1
    newSimplexMatrix(:,j)=simplexMatrix(:,j)-simplexMatrix(:,i);
end
for j=i+1:size(simplexMatrix,2)
    newSimplexMatrix(:,j-1)=simplexMatrix(:,j)-simplexMatrix(:,i);
end

end
