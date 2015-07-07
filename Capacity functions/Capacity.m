function c = Capacity(P,n)
%Capacity Calculates the symplectic capacity of the body given in P.
% Parameters
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix.
%   n - The half dimension of the body.


%parameters
iterations=1;
m=15; %Number of subdivisions of the [0:1]-interval
eps = 1e-7; %tolerance/ exactness
minAction=flintmax;

for itr=1:iterations
    %Following lines compute the matrix "A_2n" (see paper sec. 2.1 equation (2.5))
    %Matrix should be removed, and calculations should be done directly
    %to improve the efficiency of the program.
    tic
    mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)];
    toc

    A2n = zeros(2*m*n); %big zeros block
    %A2n = sparse(2*m*n,2*m*n); %%sparse matrix instead of big zeroes block
    for i=1:m
        for j=(i+1):m
            A2n(2*n*(i-1)+1:2*n*i,2*n*(j-1)+1:2*n*j)=mJ2n;
        end
    end
    %disp(A2n);

    %Generate a random initial path which fulfils the constraints
    l=0;
    while (abs(l)<eps) %look at last line in this 'while'-loop...
        x0=2*(0.5-rand(2*m*n,1)); %random numbers between -1 and 1
        s=zeros(2*n,1);
        for i=1:2*n
            s(i)=sum(x0(i:2*n:2*m*n)); %summing the i-th coordiantes, i=1,..,2n
            for j=0:m-1
                x0(j*2*n+i)=x0(j*2*n+i)-s(i)/m; %now, sum of vectors is 0
            end
        end
        l =x0'*A2n*x0; %this is needed for the last constraint
    end

    if(l<0) %in this case, the path has to be inverted
        x0p=x0;
        for i=0:m-1
            j=m-1-i;
            x0(2*n*i+1:2*n*(i+1))=x0p(2*n*j+1:2*n*(j+1));
        end
    end

    %now that l is positive, rescale x such that the last constraint is
    %fulfilled.
    l =x0'*A2n*x0;
    x0=x0*m/sqrt(l);
    %disp(x0'*A2n*x0/m^2 - 1);

  %%% USING A REGULAR POLYGON IS INITIAL POINT
% % %     rad=1;
% % %     N=m-1;
% % %     theta = [0 : 2*pi/N : 2*pi];
% % %     Pos = rad * exp(1i*theta);
% % %     X = real(Pos);
% % %     Y = imag(Pos);
% % %     regpoly = [X ; Y]';
% % %     x0=reshape(regpoly',1,2*n*(N+1));
% % %     x0=x0';
% % %     %now that l is positive, rescale x such that the last constraint is
% % %     %fulfilled.
% % %     
% % %    
% % %     l =x0'*A2n*x0;
% % %     if(l<0) %in this case, the path has to be inverted
% % %         x0p=x0;
% % %         for i=0:m-1
% % %             j=m-1-i;
% % %             x0(2*n*i+1:2*n*(i+1))=x0p(2*n*j+1:2*n*(j+1));
% % %         end
% % %     end
% % %     l =x0'*A2n*x0;
% % %     x0=x0*m/sqrt(l);
% % %     %disp(x0'*A2n*x0/m^2 - 1);
% % % 
% % %     vectB=reshape(x0,[2,m*n]);
% % %     plot(vectB(1,:),vectB(2,:));
% % %     hold on
% % %     
% % %     char = ReconstructCharacteristic(x0,P,m,n);
% % %     vect=reshape(char,[2,m*n]);
% % %     plot(vect(1,:),vect(2,:));


    cond = @(x) Constraints(x,m,n,A2n);
    pf = @(x) FuncToMinimize(x,P,m,n);
    options = optimoptions('fmincon','GradObj','on','GradConstr','on');
    %options.Display = 'iter';
    options.Algorithm = 'active-set'; %% should try which works best, Maybe this is better that sqp?
    options.MaxIter=10500;
    options.TolCon= 1e-9;
    options.TolFun= 1e-9;
    options.TolX= 1e-9;
    x=fmincon(pf,x0,[],[],repmat(eye(2*n),1,m),zeros(2*n,1),[],[],cond,options);
    %options.MaxIter=1000;
    %options.Algorithm = 'interior-point';
    %x9=fmincon(pf,x,[],[],repmat(eye(2*n),1,m),zeros(2*n,1),[],[],cond,options);


    %disp('Initial value:');
    %disp(x0);
    %disp('Final value:');
    %disp(x);
    %%disp('Iterations:');
    %%disp(k);
    
%     vectB=reshape(x,[2,m*n]);
%     plot(vectB(1,:),vectB(2,:));
%     hold on
    
    char = ReconstructCharacteristic(x,P,m,n);
    vect=reshape(char,[2,m*n]);
    plot(vect(1,:),vect(2,:));
    hold on
    
    action = 2*F(x,P,m,n);
    sprintf('Minimal action for iteration %d: %8f', itr, action)
    if (action < minAction)
        minAction = action;
    end
    
    c = minAction;
end