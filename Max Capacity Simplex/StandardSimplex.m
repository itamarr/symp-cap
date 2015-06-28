%This script calculates the capacity for the standard simplex in R^4

%for i=1:10

% this next constant c is used to translate the simplex in space such that
% 0 is contained in the simplex. This is required for the Capacity
% function.
%c = rand(1)/4;
std = [zeros(1,4);eye(4)]-1/8;% - c;

% calculate the volume of the simplex and normalize by it to have volume 1.
[R,v] = convhulln(std);
std = v ^(-1/4) * std;
avg = std'*ones(5,1) / 5;
std = std - repmat(avg',5,1);
Capacity(std,2)
% calling Capacity in this way suppresses its display
%[T, cap] = evalc('Capacity(reg,2);');
%disp('shift, capcity');
%disp([c, cap]);
%end