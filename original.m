%Code the original computer program from the paper 'Numerical Analysis in
%Hamiltonian Dynamics' by Anja Going-Jaeschke.


%parameters
tic
n=2; %(1/2)-times dimension of the space
m=8; %Number of subdivisions of the [0:1]-interval
eps = 1e-7; %tolerance/ exactness
P = [zeros(1,4);eye(4)] - 1/8;


% %initial path x0 in M_m (see paper sec. 2.2, "starting point")
% %in this case it is a square path
% e1=[1;0];
% e2=[0;1];
% x0=[-e2;-e2;-e2;-e2;-e1;-e1;-e1;-e1;e2;e2;e2;e2;e1;e1;e1;e1]/sqrt(0.125); 


%Following lines compute the matrix "A_2n" (see paper sec. 2.1 equation (2.5))
%Matrix should be removed, and calculations should be done directly
%to improve the efficiency of the program.
tic
mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)];
toc

A2n = zeros(2*m*n); %big zeros block
%A2n = sparse(2*m*n,2*m*n); %%sparse matrix instead of big zeroes block
for i=1:m
    for j=(i+1):m
        A2n(2*n*(i-1)+1:2*n*i,2*n*(j-1)+1:2*n*j)=mJ2n;
    end
end
disp(A2n);

%Generate a random initial path which fulfils the constraints
l=0;
while (abs(l)<eps) %look at last line in this 'while'-loop...
    x0=2*(0.5-rand(2*m*n,1)); %random numbers between -1 and 1
    s=zeros(2*n,1);
    for i=1:2*n
        s(i)=sum(x0(i:2*n:2*m*n)); %summing the i-th coordiantes, i=1,..,2n
        for j=0:m-1
            x0(j*2*n+i)=x0(j*2*n+i)-s(i)/m; %now, sum of vectors is 0
        end
    end
    l =x0'*A2n*x0; %this is needed for the last constraint
end

if(l<0) %in this case, the path has to be inverted
    x0p=x0;
    for i=0:m-1
        j=m-1-i;
        x0(2*n*i+1:2*n*(i+1))=x0p(2*n*j+1:2*n*(j+1));
    end
end

%now that l is positive, rescale x such that the last constraint is
%fulfilled.
l =x0'*A2n*x0;
x0=x0*m/sqrt(l);
%disp(x0'*A2n*x0/m^2 - 1);

%Initial value computed. Start the minimization process!!!

k=0; %Counts the iteration 
carryOn = true;
x=x0;
while carryOn

    %Following lines of code do the calculations in "paper, sec. 2.2, part 1"
    k = k+1;
    a = (1/m^2)*(A2n+A2n')*x; %calculates ak
    y = -dF(x,m,n); %calculates yk
    aH = proj(a,m,n); %calculates akH.  "proj(*)" calls the function "proj.m"
    yp = proj(y,m,n) - (dot(aH,y)/dot(aH,aH))*aH; %yp refers to ^yk (y with a hat). The calcualtion is similar to that of akH, so the "proj(*)" function can be reused. 
    
    if(norm(yp)<eps) % ^yk is within tolerance, so practically zero. Following step 2 in section 2.2 of the paper.
        carryOn = false;
        % we are done. x solves the minimization problem.
    else
    
        
        %Following step 3 of sec. 2.2 in the paper.
        lm=sqrt(3/4)*m/sqrt(abs(yp'*A2n*yp)); %'lm' refers to lambda_max
        hlm=dot(yp,-dF(x+lm*yp,m,n)); %'hlm' refers to hk(lambda_max)
        h0=dot(yp,yp); %h0 refers to hk(0)
        if (hlm >= 0)
            l0=lm;
        else
            l0=lm*h0/(h0-hlm);
        end
        carryOn2 = true;
        kk=0;
        if mod(k,10)==0
            disp('10-th iteration');
            disp('Minimal action:');
            action=2*F(x,m,n);
            disp(action);
            disp('error:');
            disp(norm(yp));
        end
        while carryOn2

            cl0=(1/m^2)*l0^2*(yp'*A2n*yp)+1;
            xl0=x+l0*yp;
            xml0=xl0/sqrt(cl0);
            d1=F(x,m,n)-F(xl0,m,n);
            d2=F(x,m,n)-F(xml0,m,n);
            
            if(4*d2<=d1)
                l0=l0/2; %We're not done yet, so start the while loop again
                kk=kk+1;
            else
                carryOn2=false; %we're done
                x = xml0;
            end
            if (kk==1000) %This happened earlier when 'eps=1e-8', for eps=1e-5 this is not a problem
                disp('Loop is not ending!! Ending it forcefully.');
                carryOn2=false;
                carryOn=false;
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
action=2*F(x,m,n);
disp(action);
toc
