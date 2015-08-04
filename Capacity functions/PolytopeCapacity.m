%P = [0, 0; -1, 2; 0, 3; 2, 3; 2, 1];
%P = [zeros(1,4);eye(4)];% - c;
P = [0 0 0 0; 1 0 0 0; 0 1 2 0; 0 0 3 1; 0 1 1 1]
%Capacity(P,2);
%Capacity(CentralizePolytope(LagrangianProductPolytope([0,0;pi,0;0,pi],[0,0;1,0;0,1;1,1])),2

Capacity(P,2,'subintervals', 8);
%disp('now the old one');
OldCapacity(P,2,'subintervals', 64);
