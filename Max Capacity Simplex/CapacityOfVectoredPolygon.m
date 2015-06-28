function capacity = CapacityOfVectoredPolygon(X,n)
%CapacityOfVectoredPolygon Calculates the capacity of a polygon given as
%vector.
%   X is a vector such that every 2n terms determine a point in R^2n.
%   n is the half-dimension of the ambient euclidean space.

    [rows, cols] = size(X);
    if (mod(rows,2*n) ~= 0)
        disp('error in size of X, it must be a multiple of 2n');
        return
    end
    
    P = reshape(X, [2*n, rows/(2*n)])';
    capacity = Capacity(P,n);
end