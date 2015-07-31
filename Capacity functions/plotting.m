%Parameters
%n=1
xmin = -10;
xmax = 10;
ymin =-10;
ymax=10;
xres=50;
yres=50;
polarBody=[0,0;7,1;3,5];

%
XX=linspace(xmin,xmax,xres);
YY=linspace(ymin,ymax,yres);

%
barycenter=sum(polarBody,1)/3;
polarBody=polarBody-repmat(barycenter,3,1);
%g =@(x,y) G([x;y],polarBody,1);

figure
title('the triangle');
plot([polarBody(:,1);polarBody(1,1)],[polarBody(:,2);polarBody(1,2)]);
axis equal;

figure
title('function G');
Gvals=zeros(yres,xres);
for i=1:xres
    for j=1:yres
        pp=[XX(i);YY(j)];
        Gvals(j,i)=G(pp,polarBody,1);
    end
end
surf(XX,YY,Gvals);
xlabel 'x-axis';
ylabel 'y-axis';
zlabel 'values of G';
grid on;


figure
title('function dG');
dGvals=zeros(xres*yres,4);
k=1;
for i=1:xres
    for j=1:yres
        ppp=[XX(i);YY(j)];
        dGvals(k,:)=[XX(i),YY(j),dG(ppp,polarBody,1)'];
        k=k+1;
    end
end
quiver(dGvals(:,1),dGvals(:,2),dGvals(:,3),dGvals(:,4));
