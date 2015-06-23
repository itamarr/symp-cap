% 
% generate a 5-pointed polytope in R^4
% find the convex hull and its volume
% discard this polytope if its degenerate
% rescale to get volume 1
% calculate the capcity
% find maximum 
n=2;

startingSimplex = GenerateSimplex(n);
[rows, cols] = size(startingSimplex);
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