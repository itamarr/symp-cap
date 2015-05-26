% h = waitbar(0,'please wait...');
% steps =1000;
% max =0;
% err=0;
% for ii=1:50
%     run original.m;
%     if(kk==1000)
%         err=err+1;
%         n=norm(yp);
%         if max<n
%             max=n;
%         end
%     end
%     if(mod(steps,5)==0)
%         waitbar(ii/steps);
%     end
% end
% max
% err


Data=zeros(29,4);
%for z=1:9
    n=2;
    m=160;
    
for p=2:30
    run original.m;
    analytic = ((2*gamma(2-1/p))^(2*n))/gamma(1+2*n*(1-1/p)); %Volume of the L-q-norm Ball
    %analytic = 4*gamma(2-1/p)^2/(gamma(3-2/p));
    err = abs((analytic-action)/analytic); %this only has a meaning for n=1, where the volume and the capacity coincide
    Data(p-1,:)=[p,analytic,action,err];
end
%end

figure
plot(Data(:,1),Data(:,2)) %this is the volume
hold on
plot(Data(:,1),Data(:,3),'*') %these are the capacities
