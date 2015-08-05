%
% The simplex we ran on is the convex hull of the next vertices represented
% as columns.
% These are used 
%
% [ 1,    0, -a2_3/a1_2,                                 -a2_4/a1_2]
% [ 0, a1_2,       a1_3,                                       a1_4]
% [ 0,    0,          1,                                          0]
% [ 0,    0,          0, a3_4 - (a1_3*a2_4)/a1_2 + (a1_4*a2_3)/a1_2]


results1 = dlmread('Overnight\20150804 1040 1\Simplex_Capacities_Normalized.csv');
results1 = results1(2:end,:);

vect = -5:0.5:5;
figure
title('parameters a1_2, a1_3')
bar3(vect, results1);

results2 = dlmread('Overnight\20150804 1040 2\Simplex_Capacities_Normalized.csv');
results2 = results2(2:end,:);


figure
bar3(vect, results2);

vect = -5:0.5:3;
results3 = dlmread('Overnight\20150804 1040 3\Simplex_Capacities_Normalized.csv');
results3 = results3(2:end,:);


figure
bar3(vect, results3);
