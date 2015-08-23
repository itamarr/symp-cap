function [c, char, udot] = Capacity(P,n, varargin)
%Capacity Calculates the symplectic capacity of the body given in P.
% Parameters
%   P - The convex hull of the body whose capacity we wish to calculate. 
%   It should have the form of a k-2*n matrix where k is the number of vertices
%   in the convex hull.
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
%   'MaxSubintervals' - The algorithm increases the number of subintervals
%   until it reaches this value. This defaults to 400.
%   'MinimumSubintervals' - The algorithm increases its precision after
%   this number is reached. This defaults to 30.
%   'ExponentBound' - The algorithm increases the number of subintervals by
%   2 in each step. This parameter bounds the number of such steps. This
%   defaults to 4.
%   'MinkSum' - When this parameter is turned on, the algorithm will use a
%   different set of G, and F functions. WARNING: This is still work in
%   progress and is known to be wrong. This defaults to 'off'.
%   'Epsilon' - The radius of the ball to add to the body for the Minkowski
%   sum options. This defaults to 0.25.
%   'DisplayIter' - Turn this on to display the results of each iteration
%   of fmincon. This defaults to 'off'.
%   'TolIter' - This tolerance value helps the algorithm decide to stop
%   iterating fmincon with finer subdivision of intervals. If two
%   consecutive calls to fmincon is <= this parameter, the algorithm stops.
%   This defaults to 8e-4.

tic

udot = 0;
char = 0;

% This next part sets the default parameters of the function and reads
% whatever options were input to the function.
funcParameters = struct('subintervals', 4, 'maxsubintervals', 400, 'minimumsubintervals', 30, 'exponentbound', 4, 'plotchar', 'off', 'plotudot', 'off', 'iterations', 1, 'minksum', 'off', 'epsilon', 0.25, 'startingtraj', 0, 'toliter', 0.8*1e-3, 'displayiter', 'off');
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

%parameters
iterations = funcParameters.iterations;
m = funcParameters.subintervals; %Number of subdivisions of the [0:1]-interval
eps = 1e-7; %tolerance/ exactness
minAction=flintmax;

%Following lines compute the matrix "A_2n" (see paper sec. 2.1 equation (2.5))
%Matrix should be removed, and calculations should be done directly
%to improve the efficiency of the program.
mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)];

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
    if (strcmpi(funcParameters.displayiter, 'on'))
        options.Display = 'iter';
    end
    options.Algorithm = 'active-set'; %% should try which works best, Maybe this is better that sqp?
    options.MaxIter = 8000;
    options.TolCon = 1e-9;
    options.TolFun = 1e-3;
    options.TolX = 1e-3;

    % fmincon is Matlab's function that finds minimum under constraints.
    [prevx, prevval] = fmincon(pf,x0,[],[],repmat(eye(2*n),1,m),zeros(2*n,1),[],[],cond,options);
    TolIter = funcParameters.toliter;
    for k = 1 : funcParameters.exponentbound
        tmpprevx = SubdividePath(prevx,n);
        prevx = tmpprevx;
        m = 2*m;
        disp(m);
        %%%generate A2n for new value of m
        A2n = zeros(2*m*n); %big zeros block
        for i = 1:m
            for j = (i + 1):m
                A2n(2*n*(i-1)+1:2*n*i,2*n*(j-1)+1:2*n*j)=mJ2n;
            end
        end
        % increase the precision of the algorithm of m is big enough or
        % we've had enough iterations.
        if (k >= 2 || m >= funcParameters.minimumsubintervals)
            options.TolFun = 1e-6;
            options.TolX = 1e-6;
        end
        
        %%%  reassign functins with tne new value of m

         cond = @(x) Constraints(x,m,n,A2n);
         if (strcmpi(funcParameters.minksum, 'off'))
             pf = @(x) FuncToMinimize(x,P,m,n);
         else
             pf = @(x) FuncToMinimizeMinkSum(x,P,m,n, funcParameters.epsilon);
         end
                
        [newx, newval] = fmincon(pf,prevx,[],[],repmat(eye(2*n),1,m),zeros(2*n,1),[],[],cond,options);
        if abs(newval - prevval) <= TolIter
            prevx=newx;
            break
        else
            prevval = newval;
            l = newx'*A2n*newx;
            newx = newx*m/sqrt(l);
            prevx = newx;
        end
        
    end
    udot = prevx;

    if (strcmpi(funcParameters.plotchar, 'on'))
         char = ReconstructCharacteristic(udot,P,m,n);
         vect = reshape(char,[2,m*n]);
         scatter(vect(1,:),vect(2,:));
         hold on
    end

    if (strcmpi(funcParameters.plotudot, 'on'))
         vectU = reshape(udot, [2,m*n]);
         scatter(vectU(1,:), vectU(2,:));
         hold on
    end

    action = 2*F(udot,P,m,n);
    sprintf('Minimal action for iteration %d: %8f', itr, action)
    if (action < minAction)
        minAction = action;
    end
    
    c = minAction;
    
    toc
end
