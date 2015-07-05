function result = GenerateSimplexParametrization(n)
% GenerateSimplexParamertization Generates a parametrization for the space
% of simplexes in R^{2n}, in symbolic form. The ai_j symbols represent the
% symplectic area of 2-dim faces of the simplex. 
% To use the result, call subs(result, 'ai_j', <value>), where instead of 
% i and j you should place the actual indices, and <value> should be a real number.
%   n - half the dimension of the simplex.
%   result - Here will be returned a 2n-by-2n matrix, whose terms are
%   symbols of real valued parameters which are the symplectic areas of
%   2-dim faces of the simplex.

    A = sym('a', [2*n 2*n], 'real');
    % make A an upper triangular matrix, since we know that the i,j-th
    % terms where j>i are 0.
    A = triu(A);
    % The first couple of vectors are predetermined.
    A(1,1) = 1;
    A(2,2)= A(1,2);
    A(1,2)=0;
    
    % Generate the complex structure. NOTE: Here the coordinates are given
    % by pairs q_i,p_i (instead of q_1,...,q_n,p_1,...,p_n).
    % MAKE SURE THIS COINCIDES WITH HOW WE COMPUTE CAPACITIES
    J1 = [0 1;-1 0];
    Jn = kron(eye(n), J1);
    G = abs(Jn);
    % this determines the order in which we set the coordinates of the
    % current vectors
    indices = G * [1:2*n]';

    for j = 3:2:2*n
        % To find u_j we need to calculate \omega(i col of A, u_j), which is
        % exactly dot((i col of A), Jn * u_j)
        
        % use temporary variables for the currently calculated vectors.
        % This is needed since A is used both as the original parameters
        % matrix and as an iterative result.
        v = A(:,j);
        v2 = A(:,j+1);
        prod = Jn * A(:,j);
        prod2 = Jn * A(:,j+1);
        % param and param2 are vectors of parameters. We keep them
        % seperately because otherwise we'd run over them as we set values
        % back into A.
        param = A(:,j);
        param2 = A(:,j+1);

        for i = 1:j-1
            v(indices(i)) = solve(dot(A(:,i), prod) == A(i,j), param(indices(i)));
            v2(indices(i)) = solve(dot(A(:,i), prod2) == A(i,j+1), param2(indices(i)));
        end

        % To determine the last couple of terms in the j-th and (j+1)-st
        % vectors, we have one formula to solve:
        % a_{j,j+1} = \omega(v_j,v_{j+1})=\omega(u_j,u_{j+1})
        
        A(:,j)=v;
        A(:,j+1)=v2;
        A(j,j)=1;
        temp = solve(dot(A(:,j),Jn * A(:,j+1))==A(j,j+1),param2(j+1));
        v2(j)=0;
        v2(j+1)= temp;
        A(:,j+1)=v2;
    end

    result = A;
end
