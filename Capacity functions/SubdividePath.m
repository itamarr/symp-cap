function [ y ] = SubdividePath( tx, n )
%UNTITLED2 Takes a path: a vector consisting of m concatenated vectors in
%R^2n and returns a path of length 2m by insering the average of every two
%points between them
%   Detailed explanation goes here
x=tx';
m = uint32(size(x,2)/(2*n));
if ~isinteger(m)
    disp('length of x is not a multiple of 2n.');
end
y=zeros(1,4*n*m);
A = zeros(2,size(x,2));
A(1,:) = x;
for i=1:m-1
    A(2,2*n*(i-1)+1:2*n*i) = 0.5*(x(2*n*(i-1)+1:2*n*i)+x(2*n*i+1:2*n*(i+1)));
end
A(2,2*n*(m-1)+1:2*n*m) = 0.5*(x(2*n*(m-1)+1:2*n*m)+x(1:2*n));

sizesArray = repmat(2*n,1,m);
y=cell2mat(reshape(mat2cell(A,[1,1],sizesArray),1,2*m));
y=y';
end

