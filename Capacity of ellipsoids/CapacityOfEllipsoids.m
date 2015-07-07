function c = CapacityOfEllipsoids(R,n)
%Capacity Calculates the symplectic capacity of the ellipsoid whose
%spectrum is given in R. The purpose of this is a sanity check to see that
%we get the correct analytic answer. This is done in a similar way to the 
%general case, but the h_K used is tailored for ellipsoids.
% Parameters
%   R - The spectrum of the ellipsoid in the form of an n-dimensional
%   vector in an ascending order.
%   n - The half dimension of the body.


%parameters
iterations=1;
m=12; %Number of subdivisions of the [0:1]-interval
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





    cond = @(x) Constraints(x,m,n,A2n);
    pf = @(x) FuncToMinimizeEllipsoids(x,R,m,n);
    options = optimoptions('fmincon','GradObj','on','GradConstr','on');
    %options.Display = 'iter';
    options.Algorithm = 'active-set'; %% should try which works best, Maybe this is better that sqp?
    options.MaxIter=1000;
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
    action = 2*F_Ellipsoids(x,R,m,n);
    sprintf('Minimal action for iteration %d: %8f', itr, action)
    if (action < minAction)
        minAction = action;
    end
    
    c = minAction;
end