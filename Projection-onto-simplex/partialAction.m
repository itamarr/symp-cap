function [ action ] = partialAction(char,a,b,n)
% This function computes the action of a subinterval of a characteristic.
% Since the characteristic vector "char" gives us points, we compute the
% tangent vectors by taking a difference, and the points by taking an
% average. 

action = 0;
J = [zeros(n),eye(n);-eye(n),zeros(n)];

for i=a:b
    u = J*(char(i*2*n+1:(i+1)*2*n)-char((i-1)*2*n+1:i*2*n))';
    v = char(i*2*n+1:(i+1)*2*n)+char((i-1)*2*n+1:i*2*n);
    v = 0.5*v;
    action = action + dot(-u',v);
    i = i+1;
end
action = 0.5*action;
end
