overnightUdots = dlmread('Results/20150804 1040 2/udots.csv');

% charNum is the row number to read from the file.
% Each row stands for a new characteristic.
charNum = 100;

n=2;
params = GenerateSimplexParametrization(n);

% set all the parameters to 1, beside your favorite 2 parameters.
for i = 1:2*n
    for j = i+1:2*n
        if (~(i == 1 && j == 2) && ~(i == 1 && j == 4))
            params = subs(params, sprintf('a%d_%d', i,j), 1);
        end
    end
end

a = overnightUdots(charNum,1);
b = overnightUdots(charNum,2);

P = zeros(2*n + 1, 2*n);
temp = subs(params, 'a1_2',a);
temp = subs(temp, 'a1_4', b);
P(1:end-1,:) = temp;

udot = double(overnightUdots(charNum, 5:end)');

m = real(size(udot, 1) / (2*n));

char = ReconstructCharacteristic(udot, P, m, n);

[faceArray,alpha] = prSimplex(2,P',char,0.01)
