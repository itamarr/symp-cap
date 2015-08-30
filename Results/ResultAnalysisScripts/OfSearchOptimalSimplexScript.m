% This script goes over the results of the script SearchOptimalSimplex.m
% and tries to reproduce the abnormal results.
n = 2;
numOfReprod = 30; % the number of times we want to reproduce the results
viterboConstPlusError = sqrt(2) - 0.05;

% Get the list of values who refute Viterbo's conjecture
abnormalNormCapResults = volNormCapacityResults(find(volNormCapacityResults > viterboConstPlusError));
abnormalSimplices = simplices(find(volNormCapacityResults > viterboConstPlusError));

redoCapacities = zeros(numOfReprod, size(abnormalNormCapResults, 2));
redoNormResults = zeros(numOfReprod, size(abnormalNormCapResults, 2));

h = waitbar(0, 'initializing waitbar...');
numOfOddResults = size(abnormalNormCapResults,2);

% Plot the initial odd results, together with the thershold determined by
% Viterbo's conjecture.
figure
hold on
scatter(1:numOfOddResults, abnormalNormCapResults);
plot(1:numOfOddResults, repmat(sqrt(2),1, numOfOddResults));

for j = 1:numOfReprod
    for i = 1:numOfOddResults
        P = cell2mat(abnormalSimplices(i));
        [A, vol] = convhulln(P);

        waitbar((j*numOfOddResults + i) / (numOfOddResults * numOfReprod), h, sprintf('Plot no. %d. Calculating capacity for simplex no. %d', j, i));
        [c, char, udot] = Capacity(P, n);

        redoCapacities(j,i) = c;
        redoNormResults(j,i) = c/ (vol^(1/n));
    end
    
    % Plot the reproduced results.
    scatter(1:numOfOddResults, redoNormResults(j,:));
end
close(h)

averages = sum(redoNormResults) / numOfReprod;
find(averages > viterboConstPlusError)
