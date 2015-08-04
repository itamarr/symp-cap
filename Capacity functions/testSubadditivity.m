tic
% standardSimplex=subdividebyhyperplane(0,0,1,1,0,0);
% maxSymmSimplex=subdividebyhyperplane(1,1,1,1,1,1);
%time in h:m:s
h=0;
m=10;
s=0;
t=h*3600+m*60+s;
%randSimplices=cell(0);
while toc<t
    A=10*rand;
    B=10*rand;
    C=10*rand;
    D=10*rand;
    E=10*rand;
    F=10*rand;
    try
        randSimplices{end+1}=subdividebyhyperplane(A,B,C,D,E,F);
    catch E
        warning('horrible stuff happened');
    end
end
