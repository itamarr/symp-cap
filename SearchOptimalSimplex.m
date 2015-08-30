eps = 1e-6;
n = 2;
iterations = 380;
maxEntryValue = 100;
currentIteration = 1;

% result vectors
capacityResults = zeros(1, iterations);
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
    
    % Calculate the capacity of the simplex, safely.
    try
        waitbar(currentIteration / iterations, h, sprintf('Calculating capacity for simplex no. %d', currentIteration));
        [c,char, udot] = Capacity(P,n);
    catch ME
        warning('Encountered warning while calculating capacity. May god help us all: %s\n', ME.msgtext);
        fprintf(logFile, 'Encountered warning while calculating capacity. May god help us all: %s\n', ME.msgtext);
    end
    
    % Save the results to file, including the ratio between the capacity
    % and the n-th root of the simplex's volume.
    capacityResults(currentIteration) = c;
    volNormCapacityResults(currentIteration) = c/ (vol^(1/n));
    udots(currentIteration) = mat2cell(udot,size(udot,1));
    
    dlmwrite('Random_Simplex_Capacities_Results.csv', [c, volNormCapacityResults(currentIteration), udot'], '-append');
    
    
    currentIteration = currentIteration + 1;
end

% get the list of values who refute Viterbo's conjecture
% abnormalNormCapResults = volNormCapacityResults(find(volNormCapacityResults > sqrt(2)-0.05))


% dispose of trash and save workspace
close(h);
fclose(logFile);
save('Random_Simplex_Capacities_Results_Workspace');

