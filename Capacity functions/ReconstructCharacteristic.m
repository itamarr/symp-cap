function [ Char ] = ReconstructCharacteristic( x, P, m, n )
%RECONSTRUCTCHARACTERISTIC Summary of this function goes here
%   Detailed explanation goes here

%implementation of the reconstruction from the paper
% note that in the paper the elements of (R^n)^m are indexed as x_0 ...
% x_{m-1}
% A B C stand for each term in formula (2.14)
lambda = F(x,P,m,n);
mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)]; %'mJ2n' is -J_2n.

matX = reshape(x,[2*n,m]);

v = zeros(size(x));



sumOfGrads = 0;
for j=0:m-1
    sumOfGrads = sumOfGrads + dG(mJ2n*matX(:,j+1),P,n);
end
B = (1/(sqrt(lambda)*m))*sumOfGrads;

tmpV=zeros(2*n,1);
for i=1:m-1
    for j=0:i-1
        tmpV=tmpV+matX(:,j+1);
    end
end
C = (-2*sqrt(lambda)/(m^2))*tmpV;

v(1:2*n) = B + C;

for k=1:m-1
    A = (2*sqrt(lambda)/m)*sum(matX(:,1:k),2); 
    v(k*2*n+1:(k+1)*2*n) = A+B+C;
end



Char = reshape(v,[1,2*n*m]);
end

