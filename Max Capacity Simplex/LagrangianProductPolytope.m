function [ extLagP ] = LagrangianProductPolytope( P1, P2 )
%LAGRANGIANPRODUCT gets two n-dim Polytopes and returns the extreme points
%of the simplex resulting in their lagrangian product (in coordinates
%q1,..,qn,p1,...,pn
%   Detailed explanation goes here
if (size(P1,2) ~= size(P2,2))
    disp ('Simplices must be of the same dimension');
end

lagProd = zeros(size(P1,1)*size(P2,1),size(P1,2)+size(P2,2)); %get [row1*row2, n] empty matrix
for i=1:size(P1,1)
    for j=1:size(P2,1)
        lagProd(((i-1)*size(P2,1)+j),:)=[P1(i,:),P2(j,:)]; %concatenate the appropriate vectors
    end
end

extLagP = ExtremePointsOnly(lagProd);

end

