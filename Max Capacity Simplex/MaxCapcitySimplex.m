% This script finds the 5-pointed polytope in R^4 of volume 1 with maximal
% capacity. This is done by generating a random non-degenerate polytope,
% then using fmincon on minus the capacity function, with the constraint
% that the volume of each step is 1.

n=2;

startingSimplex = GenerateSimplex(n);
[rows, cols] = size(startingSimplex);

% break the starting simplex into one big vector since this is the input
% for fmincon
x0 = reshape(startingSimplex, [rows*cols, 1]);
conditions = @(x) VolumeConstraint(x,n);
capacityFun = @(x) -1 * CapacityOfVectoredPolygon(x,n);
options = optimoptions('fmincon');
%options.Display = 'iter';
options.Algorithm = 'active-set'; %% should try which works best, Maybe this is better that sqp?
options.MaxIter=10500;
%options.TolCon= 1e-7;
%options.TolFun= 1e-8;
options.TolX= 1e-9;
x=fmincon(capacityFun,x0,[],[],[],[],[],[],conditions,options);

disp(x)

Capacity(x, n)