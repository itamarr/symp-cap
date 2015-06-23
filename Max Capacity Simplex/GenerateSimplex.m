function P = GenerateSimplex(n)
    N = 2*n;
    P = rand(N+1,N);
    % we use convhulln to compute the volumn
    [K, v] = convhulln(P);

    % make sure we don't get a degenerate polytope
    while (rank(K) < N)
        P = rand(N+1,N);
        [K, v] = convhulln(P);
    end

    % rescale the points to get volume 1
    P = v^(-1/N)*P;

    % translate the body so that it contains 0 (we need this for the
    % calculations later)
    avg = P'*ones(N+1,1) / (N+1);
    P = P - repmat(avg',N+1,1);

end