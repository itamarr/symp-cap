function result = GenerateSimplexParamertization(n)
    t = sym('t', [2*n 2*n], 'real');
    t = triu(t);
    t(1,1) = 1;
    t(2,2)= t(1,2);
    t(1,2)=0;
    
    J1 = [0 1;-1 0];
    Jn = kron(eye(n), J1);
    G = abs(Jn);
    indices = G * [1:2*n]';

    v = t(:,3);
    rep = Jn * t(:,3);
    u = t(:,3);
    % To find u_j we need to calculate \omega(i col of t, u_j), which is
    % exactly dot((i col of t), Jn * u_j)
    v(indices(1))=solve(dot(t(:,1), rep)==t(1,3), u(indices(1)));
    v(indices(2))=solve(dot(t(:,2), rep)==t(2,3), u(indices(2)));
    
    t(:,3)=v;
    
    v = t(:,4);
    rep = Jn * t(:,4);
    u = t(:,4);
    v(indices(1))=solve(dot(t(:,1),rep)==t(1,4), u(indices(1)));
    v(indices(2))=solve(dot(t(:,2),rep)==t(2,4), u(indices(2)));
    
    % in this example, we've gone up to i=j-2, now we need to settle the
    % rest of the matrix which is a 2-by-2 block with 2 degrees of freedon
    
    t(3,3)=1;
    t(:,4)=v;
    temp = solve(dot(t(:,3),Jn * t(:,4))==t(3,4),u(4));
    v(3)=0;
    v(4)= temp;
    t(:,4)=v;
    
    % for j=3:2:2*n
    % v = Jn * t(:,j);
    % v2 = Jn * t(:,j+1);
    % u = t(:,j);
    % rep = Jn * t(:,j);
    % for i=1:j-2
    %   v(indices(i))=solve(dot(t(:,i),rep)==t(i,j), u(indices(i)))
    %   v2(indices(i))=solve(dot(t(:,i),rep)==t(i,j+1), u(indices(i)))
    % end
    % solve(dot(t(:,j),t(:,j+1))==t(j,j+1),u(indices(j+1))
    % v(j)=0;
    % end
    
    %% Note %%
    % to substitute real values into the result, use: subs(t, 'ti_j',
    % number). The output of this command is t where the term ti_j is
    % replaced by number.
    result = t;
end
