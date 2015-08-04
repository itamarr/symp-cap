% Most general simplex: P = [1,0,0,0; 0,a,1,0; 0,b,0,d; 0,c,f/e,e; 0,0,0,0]
eps = 1e-6;
n = 2;
% Create a simplex parametrization of the correct dimensions
params = GenerateSimplexParametrization(n);

% set all the parameters to 1, beside your favorite 2 parameters.
for i = 1:2*n
    for j = i+1:2*n
        if (~(i == 1 && j == 2) && ~(i == 1 && j == 4))
            params = subs(params, sprintf('a%d_%d', i,j), 1);
        end
    end
end

% these are the a and b which we're going to substitute for parameters
aVect = -5:0.5:5;
bVect = -5:0.5:5;

% result vectors
capacityresults = zeros(size(aVect,2),size(bVect,2));
volnormcapacityresults = zeros(size(capacityresults));
udotsVect = zeros(size(capacityresults));

dlmwrite('Simplex_Capacities_Raw.csv', aVect, '-append');
dlmwrite('Simplex_Capacities_Normalized.csv', bVect, '-append');

logFile = fopen('SimplexCapacityScript.log', 'w');

i=1;
for a = aVect
    j=1;
    for b = bVect
        %P = [1,0,0,0; 0,1,a,0; 0,0,b,1; 0,1,1,1; 0,0,0,0] 
        
        % dirty fix - we arbitrarily decided that a1_2 is not zero... and
        % sometimes we want it to be zero. so.. emr... never mind and
        % ignore it.
        if (a == 0)
            warning('we hit a parameter that gives zero volume simplex for a=%d, b=%d; ignore and continue.', a, b);
            fprintf(logFile, 'we hit a parameter that gives zero volume simplex for a=%d, b=%d; ignore and continue.', a, b);
            j = j+1;
            continue;
        end
        P = zeros(2*n + 1, 2*n);
        temp = subs(params, 'a1_2',a);
        temp = subs(temp, 'a1_4', b);
        P(1:end-1,:) = temp;
        %n = size(P,2)/2 % no need to calculate n, it's fixed at the
        %beginning
        vol = abs(det(P(1:end-1,:))/factorial(2*n))
        if (vol < eps)
            warning('skipping simplex with volume close to zero, with parameters a=%f, b=%f', a, b);
            fprintf(logFile, 'skipping simplex with volume close to zero, with parameters a=%f, b=%f', a, b);
            j = j + 1;
            continue
        end

        try
            [c,char, udot] = Capacity(P,n);
        catch ME
            warning('Encountered warning while calculating capacity. May god help us all: %s', ME.msgtext);
            fprintf(logFile, 'Encountered warning while calculating capacity. May god help us all: %s', ME.msgtext);
        end
        udotsVect(i,j) = udotsVect(i, j);
        capacityresults(i,j) = c
        volnormcapacityresults(i,j) = c/(vol^(1/n))
        
        % this writes the udot result of Capacity function to file, for
        % future generations
        dlmwrite('udots.csv', [a,b,c,volnormcapacityresults(i,j), udot'], '-append'); 

        j = j + 1;
    end
    
    % this saves the raw data of the capacity function before and after
    % normalization
    dlmwrite('Simplex_Capacities_Raw.csv', capacityresults(i,:), '-append');
    dlmwrite('Simplex_Capacities_Normalized.csv', volnormcapacityresults(i,:), '-append');
    i = i + 1;
end

fclose(logFile);


