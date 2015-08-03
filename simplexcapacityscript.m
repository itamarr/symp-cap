capacityresults = zeros(20,20);
volnormcapacityresults = [zeros(20,20)];
% Most general simplex: P = [1,0,0,0; 0,a,1,0; 0,b,0,d; 0,c,f/e,e; 0,0,0,0]
eps = 1e-6;

i=1;
for a = -0.8:0.2:0.8
    j=1;
    for b = -0.8:0.2:0.8
        P = [1,0,0,0; 0,a,1,0; 0,b,0,1; 0,1,1,1; 0,0,0,0] 
        n = size(P,2)/2
        vol = abs(det(P(1:end-1,:))/factorial(2*n))
        if (vol < eps)
            warning('skipping simplex with volume close to zero, with parameters a=%f, b=%f', a, b);
            continue
        end

        [c,char] = Capacity(P,n);
        capacityresults(i,j) = c
        volnormcapacityresults(i,j) = c/(vol^(1/n))

        j = j + 1;
    end
    i = i + 1;
end

dlmwrite('SimplexCapacities.csv', a, '-append');
dlmwrite('SimplexCapacities.csv', b, '-append');
dlmwrite('SimplexCapacities.csv', capacityresults, '-append');
dlmwrite('SimplexCapacities.csv', volnormcapacityresults, '-append');

