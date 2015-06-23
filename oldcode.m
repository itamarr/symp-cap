
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
