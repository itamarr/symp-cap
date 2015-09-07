eps = 1e-6;
n = 2;

% Create a simplex parametrization of the correct dimensions
params = GenerateSimplexParametrization(n);

ai = 2;
aj = 3;
bi = 1;
bj = 4;

fixedParam = 1;

% set all the parameters to 1, beside your favorite 2 parameters.
for i = 1:2*n
    for j = i+1:2*n
        if (~(i == ai && j == aj) && ~(i == bi && j == bj))
            params = subs(params, sprintf('a%d_%d', i,j), fixedParam);
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
        % dirty fix - we arbitrarily decided that a1_2 is not zero... and
        % sometimes we want it to be zero. so.. emr... never mind and
        % ignore it.
        if (ai == 1 && aj == 2 && a == 0)
            warning('we hit a parameter that gives zero volume simplex for a=%d, b=%d; ignore and continue.\n', a, b);
            fprintf(logFile, 'we hit a parameter that gives zero volume simplex for a=%d, b=%d; ignore and continue.\n', a, b);
            j = j+1;
            continue;
        end
        P = zeros(2*n + 1, 2*n);
        temp = subs(params, sprintf('a%d_%d', ai, aj), a);
        temp = subs(temp, sprintf('a%d_%d', bi, bj), b);
        P(1:end-1,:) = temp;

        vol = abs(det(P(1:end-1,:))/factorial(2*n))
        if (vol < eps)
            s = sprintf('skipping simplex with volume close to zero, with parameters a%d_%d=%f, b%d_%d=%f\n', ai, aj, a, bi, bj, b);
            warning(s);
            fprintf(logFile, s);
            j = j + 1;
            continue
        end
        
        % normalize to volume 1
        P = P.* vol^(-1/(2*n));

        try
            [c,char, udot] = Capacity(P,n);
        catch ME
            warning('Encountered warning while calculating capacity. May god help us all: %s\n', ME.message);
            fprintf(logFile, 'Encountered warning while calculating capacity. May god help us all: %s\n', ME.message);
        end
        udotsVect(i,j) = udotsVect(i, j);
        capacityresults(i,j) = c

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


