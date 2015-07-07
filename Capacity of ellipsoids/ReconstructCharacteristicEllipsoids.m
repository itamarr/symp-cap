function [ Char ] = ReconstructCharacteristicEllipsoids( x, R, m, n )
%RECONSTRUCTCHARACTERISTIC Summary of this function goes here
%   Detailed explanation goes here
lambda = F_Ellipsoids(x,R,m,n);
lambdaSquareRoot = sqrt(lambda);
mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)]; %'mJ2n' is -J_2n.

matX = reshape(x,[2*n,m]);

v = zeros(size(x));



sumOfGrads = 0;
for j=0:m-1
    sumOfGrads = sumOfGrads + dG_Ellipsoids(mJ2n*matX(:,j+1),R,n);
end
B = (1/(lambdaSquareRoot*m))*sumOfGrads;

tmpV=zeros(2*n,1);
for i=1:m-1
    for j=0:i-1
        tmpV=tmpV+matX(:,j+1);
    end
end
C = (-2*lambdaSquareRoot/(m^2))*tmpV;

%v(1:2*n) =(2*lambdaSquareRoot/m)*matX(:,1) + B + C;

for k=0:m-1
    A = (2*lambdaSquareRoot/m)*sum(matX(:,1:k+1),2); 
    v(k*2*n+1:(k+1)*2*n) = A+B+C;
end



Char = reshape(v,[1,2*n*m]);
end

