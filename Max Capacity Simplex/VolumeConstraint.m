function [c, ceq] = VolumeConstraint(X,n)
%
       ceq = [];
       c = [];
    [rows, cols] = size(X);
    if (mod(rows,2*n) ~= 0)
        disp('error in size of X, it must be a multiple of 2n');
        return
    end
    
    P = reshape(X, [2*n, rows/(2*n)])';
    [K, v] = convhulln(P);
    ceq = v -1;
end