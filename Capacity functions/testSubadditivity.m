%time in h:m:s
h=0;
m=10;
s=0;
t=h*3600+m*60+s;
randSimplices=cell(0)
numOfRandSimplices=30;
% disp('testing candidate 1');
% i=0;
% while i<numOfTriesForCandidates
%     try
%         randSimplices{end+1}=subdividebyhyperplane(3.320312570057692,5.302492307431258,4.238524246629064,0.153229743527471,3.520270034036538,3.530366695263731);
%         i=i+1;
%     catch E
%         warning('horrible stuff happened');
%     end
% end
% disp('testing candidate 2');
% i=0;
% while i<numOfTriesForCandidates
%     try
%         randSimplices{end+1}=subdividebyhyperplane(6.503699839036440,5.312926257687318,2.629042294349858,5.058798957733863,2.314754857765802,3.679614782989666);
%         i=i+1;
%     catch E
%         warning('horrible stuff happened');
%     end 
% end
disp('testing random simplices');
i=0;
while i<numOfRandSimplices
    A=10*rand;
    B=10*rand;
    C=10*rand;
    D=10*rand;
    E=10*rand;
    F=10*rand;
    try
        randSimplices{end+1}=subdividebyhyperplane(A,B,C,D,E,F);
        i=i+1;
    catch E
        %warning('horrible stuff happened');
    end
end
save data.mat randSimplices
