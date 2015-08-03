%% Parameters

%polar body
P=makeBarycentric([0,0;3,1;5,-3]);
%range and resolution for plotting
xmin=-3;
xmax=3;
xres=100;
ymin=-3;
ymax=3;
yres=100;
%range of time (smoothing parameter)
tmin=0.1;
tmax=1;
tres=2;
numOfLevelSets = 30;
%movieSpeed = 3;

%% Finding minimum
disp('Computing...');
tic
%make grid
X=linspace(xmin,xmax,xres);
Y=linspace(ymin,ymax,yres);
T=linspace(tmin,tmax,tres);

Z=cell(tres);

%go through grid points
for k=1:tres
    Z{k}=zeros(xres,yres);
    for i=1:xres 
        for j=1:yres
            x=X(i);
            y=Y(j);
            %define function to minimize
            F = @(z) G(z,P,1)+((x-z(1))^2+(y-z(2))^2)/(2*T(k));
            guess=T(k)*dG([x;y],P,1)+[x;y];
            pMin=fminsearch(F,guess,optimset('MaxFunEvals',1000,'MaxIter',1000));
            Z{k}(i,j)=F(pMin);
            clear F;
        end
    end
end
toc

%% plotting
disp('Plotting...');
tic
figure
P=[P;P(1,:)];
plot(P(:,1),P(:,2))
title 'triangle'
axis equal
% figure('name','smoothing movie')
for k=1:tres
    figure
    M=Z{k}(:);
    LevelSets=linspace(min(M),max(M),numOfLevelSets+1);
    LevelSets=LevelSets(2:end);
    contour(X,Y,Z{k},LevelSets)
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    title(strcat('Function G_t, t = ',num2str(T(k))))
    axis equal
    %pause(movieSpeed)
end
toc
