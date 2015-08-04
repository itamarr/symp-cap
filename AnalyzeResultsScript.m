overnightUdots = dlmread('Overnight/udots.csv');

n=2;
params = GenerateSimplexParametrization(n);

% set all the parameters to 1, beside your favorite 2 parameters.
for i = 1:2*n
    for j = i+1:2*n
        if (~(i == 2 && j == 3) && ~(i == 1 && j == 4))
            params = subs(params, sprintf('a%d_%d', i,j), 1);
        end
    end
end

a = overnightUdots(1,1);
b = overnightUdots(1,2);

P = zeros(2*n + 1, 2*n);
temp = subs(params, 'a2_3',a);
temp = subs(temp, 'a1_4', b);
P(1:end-1,:) = temp;

udot = overnightUdots(1, 5:end);

m = size(udot, 2) / (2*n);

char = ReconstructCharacteristic(udot, P, m, n);

[faceArray,alpha] = prSimplex(2,P',char,0.01);
