function [c, char] = OldCapacity(P,n, varargin)
%Capacity Calculates the symplectic capacity of the body given in P.
% ParametersO
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix.
%   n - The half dimension of the body.
%   varargin - This allows for variable input arguments.
%   Usage: Capacity(P,n,'option', value,...); Where 'option' is one of the
%   available options below, and value is the desired value for the option.
%   The options are:
%   'SubIntervals' - The number of subintervals, i.e., m. This should be a
%   number. Default is 60.
%   'PlotChar' - Controls whether or not to rebuild the characteristic and
%   plot it. Values should be 'on' or 'off'. Default 'off'.
%   'PlotUDot' - Controls whether or not to plot the derivative of the
%   characteristic. Values should be 'on' or 'off'. Default 'off'.
%   'Iterations' - Controls the number of iterations for the capacity
%   function. Values should be number. Default is 1.
%   

tic

char = 0;

funcParameters = struct('subintervals', 60, 'plotchar', 'off', 'plotudot', 'off', 'iterations', 1, 'minksum', 'off', 'epsilon', 0.25, 'startingtraj', 0);
optNames = fieldnames(funcParameters);

if (round(nargin/2) ~= nargin/2)
    error('Input format is <field>, <value> pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1}); %# make case insensitive

   if any(strcmp(inpName,optNames))
      funcParameters.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

if (ischar(P))
    P = str2num(P);
end

if (ischar(n))
    n = str2double(n);
end

C = centroid(P);
P = P - repmat(C,size(P,1),1);
P=P;

%parameters
iterations = funcParameters.iterations;
m = funcParameters.subintervals; %Number of subdivisions of the [0:1]-interval
eps = 1e-7; %tolerance/ exactness
minAction=flintmax;

%Following lines compute the matrix "A_2n" (see paper sec. 2.1 equation (2.5))
%Matrix should be removed, and calculations should be done directly
%to improve the efficiency of the program.
%global mJ2n
%if (isempty(mJ2n) || size(mJ2n, 1) ~= 2*n)
    mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)];
%end

for itr=1:iterations

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
    if (funcParameters.startingtraj == 0)
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
        l = x0'*A2n*x0;
        x0 = x0*m/sqrt(l);
    else
        x0 = funcParameters.startingtraj;
    end

    cond = @(x) Constraints(x,m,n,A2n);
    if (strcmpi(funcParameters.minksum, 'off'))
        pf = @(x) FuncToMinimize(x,P,m,n);
    else
        pf = @(x) FuncToMinimizeMinkSum(x,P,m,n, funcParameters.epsilon);
    end
        
    options = optimoptions('fmincon','GradObj','on','GradConstr','on');
    %options.Display = 'iter';
    options.Algorithm = 'active-set'; %% should try which works best, Maybe this is better that sqp?
    options.MaxIter = 10000;
    options.TolCon = 1e-5;
    options.TolFun = 1e-3;
    options.TolX = 1e-3;
    x1=fmincon(pf,x0,[],[],repmat(eye(2*n),1,m),zeros(2*n,1),[],[],cond,options);
    options.TolCon = 1e-9;
    options.TolFun = 1e-6;
    options.TolX = 1e-4;
    %l =x1'*A2n*x1; %rescale to fit constraint
    %x1=x1*m/sqrt(l);
    x=fmincon(pf,x1,[],[],repmat(eye(2*n),1,m),zeros(2*n,1),[],[],cond,options);
%     options.TolFun = 1e-6;
%     options.TolX = 1e-4;
%     l =x2'*A2n*x2; %rescale to fit constraint
%     x2=x2*m/sqrt(l);
%     x=fmincon(pf,x2,[],[],repmat(eye(2*n),1,m),zeros(2*n,1),[],[],cond,options);


    if (strcmpi(funcParameters.plotchar, 'on'))
         char = ReconstructCharacteristic(x,P,m,n);
         vect=reshape(char,[2,m*n]);
         scatter(vect(1,:),vect(2,:));
         hold on
    end

    if (strcmpi(funcParameters.plotudot, 'on'))
         vectU = reshape(x, [2,m*n]);
         scatter(vectU(1,:), vectU(2,:));
         hold on
    end

    action = 2*F(x,P,m,n);
    sprintf('Minimal action for iteration %d: %8f', itr, action)
    if (action < minAction)
        minAction = action;
    end
    
    c = minAction;
    
    toc
end