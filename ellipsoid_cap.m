function [r] = ellipsoid_cap(A,c)
%[r] = ellipsoid_cap(A,c)
% inputs:
% ---------
% 2nx2n matrix A and 2nx1 vector c. Output the capacity of 
%
% outputs:
% ---------
% ellipsoid (x-c)' A (x-c) <= 1.

if size(A,1) == size(A,2),
    n = (1/2)*size(A,1);
else
    return
    display('Not even dimensional!!');    
end

if size(c,1) == 2*n,
    n = n;
else
    display('Not even dimensional!!');
    return
end

J = [zeros(n),-eye(n);eye(n),zeros(n)];

D = 2*eye(2*n);  

B = D*A;

S=eig(J*B); 

for k=1:2n
    S(k) = abs(S(k));
end

l = min(S)

r = 2*pi / l

end
