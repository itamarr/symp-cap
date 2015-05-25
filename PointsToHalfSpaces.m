function [ normals, coeff ] = PointsToHalfSpaces( P )
%PointsToHalfSpaces This function calculates the half spaces that determine
%the facets of the given polytope in the form of an equation
%normals*X=coeff
%   We calculate the convex hull of the polytope, iterate its facets,
%   calculate the normal to these facets and then calculate the
%   coefficients that determine the half space by taking the dot product of
%   the normal with the first vertex of the facet.
%
%   NOTE: We don't check for orientation issues here. We should check if
%   this is at all important.

    % convhulln returns a p by n matrix where n is the dim of the space
    % and p is the number of triangular facets of the convex hull of the
    % points in P.
    K = convhulln(P);
    normals = zeros(size(K));
    coeff = zeros(size(K,1),1);
    
    for j = 1:size(P,1)
        % this gives us the veritces of the j-th facet
        facet = P(K(j,:),:); 
        sizeOfFacetMatrix = size(facet);
        
        % we calculate linear subspaces determined by the vertices of the
        % facet. This is done by taking the difference between the first
        % vertex and every other vertex.
        % To do this efficiently we calculate this using matrices.
        D = -eye(sizeOfFacetMatrix(1));
        D(:,1) = ones(1,sizeOfFacetMatrix(1));
        parallelVectors = D*facet;
        
        % The normal to the j-th facet is then the null space of the matrix
        % whose rows are the vectors.
        normals(j,:) = (null(parallelVectors(2:end,:)))';
        
        % the dot product of the first vertex of the facet with the normal
        % to the facet gives you the coefficient that determines the half
        % plane corresponding to the facet.
        coeff(j) = dot(normals(j,:), facet(1,:));
        
        % easy fix for orientation... just invert the signs such that the
        % coefficients are positive.
        if (coeff(j) < 0)
            normals(j,:) = -1 * normals(j,:);
            coeff(j) = -1 * coeff(j);
        end
    end
end

