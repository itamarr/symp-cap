%P = [0, 0; 1, 0; 0, 1] - 1/3;
%P = [zeros(1,4);eye(4)]-1/5;% - c;
%Capacity(P,2);
Capacity(CentralizePolytope(LagrangianProductPolytope([0,0;pi,0;0,pi],[0,0;1,0;0,1;1,1])),2)
