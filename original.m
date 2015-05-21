%Code the original computer program from the paper

%parameters
tic
n=1; %(1/2)-times dimension of the space
m=16; %Number of subdivisions of the [0:1]-interval
np=2*n;
eps = 1e-8; %tolerance/ exactness
e1=[1;0];
e2=[0;1];
x0=[e1;e1;e1;e1;e2;e2;e2;e2;-e1;-e1;-e1;-e1;-e2;-e2;-e2;-e2]/sqrt(.125); %THIS CHOICE IS INVALID!! replace by a valid x0

%Matrix should be removed, and calculations should be done directly
%to improve the efficiency of the program.
mJ2n = [zeros(n),eye(n);-eye(n),zeros(n)];
A2n = zeros(2*m*n);
for i=1:m
    for j=(i+1):m
        A2n(np*(i-1)+1:np*i,np*(j-1)+1:np*j)=mJ2n;
    end
end
%disp(A2n);

k=0;
carryOn = true;
x=x0;
while carryOn
    k = k+1;
    a = (1/m^2)*(A2n+A2n')*x; %calculates ak
    y = -dF(x,m,n); %calculates yk
    aH = proj(a,m,n);
    yp = proj(y,m,n) - (dot(aH,y)/dot(aH,aH))*aH;
    if(norm(yp)<eps) % y is zero
        carryOn = false;
        % we are done. x solves the minimization problem
    else
        lm=sqrt(3/4)*m/sqrt(abs(yp'*A2n*yp));
        hlm=dot(yp,-dF(x+lm*yp,m,n));
        h0=dot(yp,yp);
        if (hlm >= 0)
            l0=lm;
        else
            l0=lm*h0/(h0-hlm);
        end
        carryOn2 = true;
        while carryOn2
            cl0=(1/m^2)*l0^2*(yp'*A2n*yp)+1;
            xl0=x+l0*yp;
            xml0=xl0/sqrt(cl0);
            d1=F(x,m,n)-F(xl0,m,n);
            d2=F(x,m,n)-F(xml0,m,n);
            if(4*d2<=d1)
                l0=l0/2; %We'e not done yet, so start the while loop again
            else
                carryOn2=false; %we're done
                x = xml0;
            end
        end
    end
end
disp('Initial value:');
disp(x0);
disp('Final value:');
disp(x);
disp('Iterations:');
disp(k);
disp('Minimal action:');
disp(2*F(x,m,n));
toc
