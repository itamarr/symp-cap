% This script calculates c(K-K) / c(K) for K the standard simplex of
% dimension 2*n, and n=2,...,50.
% The results will be saved to a *.csv file, together with the output of
% Capacity.

%% Initial parameters
maxDim = 50;
capacityResults = zeros(2, maxDim);
ratios = zeros(1, maxDim);
udots = cell(size(capacityResults));
simplices = cell(size(capacityResults));
udotK = zeros(512,1);
udotDiff = zeros(512,1);
P = zeros(10100, 100);
K = [zeros(1,100); eye(100)];

logFile = fopen('MinkowskiDiffRatio.log', 'w');

h = waitbar(0, 'initializing waitbar...');
for n=2:50
    %% Set up the polygons
    waitbar(n/maxDim, h, sprintf('Building polygon K-K for n=%d', n));
    %K = [zeros(1,2*n); eye(2*n)];
    % K-K is the convex hull of (\pm e_i, \pm (e_i-e_j))_{i,j=1}^{2n} where the e_k's
    % are the unit vectors.
    % The number of such extreme points is 2n(2n+1), so we allocate the
    % matrix for this polygon accordingly
    %P = zeros(2*n*(2*n + 1), 2*n);
    P(1:2*n,1:2*n) = eye(2*n);
    P(2*n+1:4*n,1:2*n)=-eye(2*n);

    % We calculate the rest of P in an ad hoc computation by summing the
    % i-th row with the j-th row, where i goes over the unit vectors and
    % e_j goes over the negative unit vectors.
    % The skip parameter is used to skip rows where e_i+e_j=0.
    skip = 0;
    for i=1:2*n
        for j=1:2*n
            if (i~=j)
                P(i*2*n+2*n+j-skip,:) = P(i,:)+P(2*n+j,:);
            else
                skip = skip + 1;
            end
        end
    end

    %% Calculate Capacities
    try
        waitbar(n/maxDim, h, sprintf('Calculating c(K-K)/C(K) for n=%d', n));
        [c, char, udotK] = Capacity(K(1:2*n+1, 1:2*n), n);
        capacityResults(1,n) = c;
        %udots(1, n) = mat2cell(udotK,size(udotK,1));
        fprintf(logFile, 'udotK size %d\n', size(udotK,1));
        
        [c, char, udotDiff] = Capacity(P(1:2*n*(2*n + 1), 1:2*n), n);
        capacityResults(2,n) = c;
        %udots(2, n) = mat2cell(udotDiff,size(udotDiff,1));
        fprintf(logFile, 'udotDiff size %d\n', size(udotDiff,1));
    catch ME
        warning('Encountered warning while calculating capacity. May god help us all: %s\n', ME.msgtext);
        fprintf(logFile, 'Encountered warning while calculating capacity. May god help us all: %s\n', ME.msgtext);
    end
    
    ratios(n) = capacityResults(2,n) / capacityResults(1,n);
    
    % Write the results to file. The idea of this format is that each row
    % keeps the data of a single run. Then, to differentiate between the
    % two different udots, we save its size so when we read back the data
    % we can tell how many entries did the vector have.
    dlmwrite('MinkowskiDiffRatio_Results.csv', [n,capacityResults(1,n), ...
       size(udotK,1), udotK', capacityResults(2,n), size(udotDiff,1), udotDiff'], '-append');

end

%% Clean up
% clear up the open handles and save the workspace.
fclose(logFile);
close(h);
save('MinkowskiDiffRatio_Results_Workspace');
datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')
