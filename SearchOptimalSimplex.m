eps = 1e-8;
n = 2;
iterations = 140;
maxEntryValue = 10;
capacityReruns = 5;
currentIteration = 1;

% result vectors
capacityResults = zeros(capacityReruns, iterations);
volNormCapacityResults = zeros(size(capacityResults));
udots = cell(size(capacityResults));
simplices = cell(size(capacityResults));

logFile = fopen('SimplexCapacityScript.log', 'w');

h = waitbar(0, 'initializing waitbar...');
while (currentIteration <= iterations)
    % Generate a random simplex and make sure its volume is not too small
    P = maxEntryValue.* rand(2*n+1, 2*n);
    [A, vol] = convhulln(P);
    if (vol < eps)
        s = sprintf('skipping simplex with volume close to zero\n');
        warning(s);
        fprintf(logFile, s);
        continue
    end
    
    % Save the data of the current simplex to file and to the workspace
    simplices(currentIteration) = mat2cell(P,size(P,1),size(P,2));
    dlmwrite('Random_Simplex_Capacities_Results.csv', P, '-append');
    
    for j = 1:capacityReruns
        % Calculate the capacity of the simplex, safely.
        try
            waitbar((currentIteration*capacityReruns + j) / (capacityReruns * iterations), h, sprintf('Calculating capacity for simplex no. %d, run no. %d', currentIteration, j));
            [c,char, udot] = Capacity(P,n);
        catch ME
            warning('Encountered warning while calculating capacity. May god help us all: %s\n', ME.msgtext);
            fprintf(logFile, 'Encountered warning while calculating capacity. May god help us all: %s\n', ME.msgtext);
        end

        % Save the results to file, including the ratio between the capacity
        % and the n-th root of the simplex's volume.
        capacityResults(j, currentIteration) = c;
        volNormCapacityResults(j, currentIteration) = c/ (vol^(1/n));
        udots(j, currentIteration) = mat2cell(udot,size(udot,1));
    
        dlmwrite('Random_Simplex_Capacities_Results.csv', [c, volNormCapacityResults(j, currentIteration), udot'], '-append');
    end
    
    
    currentIteration = currentIteration + 1;
end

% get the list of values who refute Viterbo's conjecture
% abnormalNormCapResults = volNormCapacityResults(find(volNormCapacityResults > sqrt(2)-0.05))


% dispose of trash and save workspace
close(h);
fclose(logFile);
save('Random_Simplex_Capacities_Results_Workspace');
datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')
