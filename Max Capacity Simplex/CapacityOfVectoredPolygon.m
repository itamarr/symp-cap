function capacity = CapacityOfVectoredPolygon(X,n)
    [rows, cols] = size(X);
    if (mod(rows,2*n) ~= 0)
        disp('error in size of X, it must be a multiple of 2n');
        return
    end
    
    P = reshape(X, [2*n, rows/(2*n)])';
    capacity = Capacity(P,n);
end