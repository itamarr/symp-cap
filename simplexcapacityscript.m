% Most general simplex: P = [1,0,0,0; 0,a,1,0; 0,b,0,d; 0,c,f/e,e; 0,0,0,0]
eps = 1e-6;

aVect = -5:0.5:5;
bVect = -5:0.5:5;

capacityresults = zeros(size(aVect,2),size(bVect,2));
volnormcapacityresults = zeros(size(capacityresults));
udotsVect = zeros(size(capacityresults));

dlmwrite('Simplex_Capacities_Raw.csv', aVect, '-append');
dlmwrite('Simplex_Capacities_Normalized.csv', bVect, '-append');

i=1;
for a = aVect
    j=1;
    for b = bVect
        P = [1,0,0,0; 0,1,a,0; 0,0,b,1; 0,1,1,1; 0,0,0,0] 
        n = size(P,2)/2
        vol = abs(det(P(1:end-1,:))/factorial(2*n))
        if (vol < eps)
            warning('skipping simplex with volume close to zero, with parameters a=%f, b=%f', a, b);
            j = j + 1;
            continue
        end

        [c,char, udot] = Capacity(P,n);
        udotsVect(i,j) = udotsVect(i, j);
        capacityresults(i,j) = c
        volnormcapacityresults(i,j) = c/(vol^(1/n))
        
        dlmwrite('udots.csv', [a,b,c,volnormcapacityresults(i,j), udot'], '-append'); 

        j = j + 1;
    end
    
    dlmwrite('Simplex_Capacities_Raw.csv', capacityresults(i,:), '-append');
    dlmwrite('Simplex_Capacities_Normalized.csv', volnormcapacityresults(i,:), '-append');
    i = i + 1;
end


