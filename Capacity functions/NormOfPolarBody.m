function [ n ] = NormOfPolarBody( P, y )
%NormOfPolarBody Calculates the norm of y in relation to the polar of P.


    %[A, b] = PointsToHalfSpaces(P);
    %%A(~any(A,2),:) = [];
    %%b(~any(b,2),:) = [];
    
    [M,I] = max(P*(y)'); %% it's faster!
    
    
%     tic
%     options = optimset('Display','none'); %% disable logging for linprog
%     maxX = linprog(-y,A,b,[],[],[],[],[],options);
%     disp('linprog');
%     toc
%     n = y * maxX;
    n = (y * (P(I(1),:)'));
end

