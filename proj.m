function [vH] = proj(v,m,n)
%PROJ Projection onto gi=0-plane
%   Projection of a vector onto the H-plane
%   described in section 2.3 of the paper, equation (2.17)

vHp=zeros(2*m*n,1);
for i=1:2*n
    s=0; %the scalar product
    for j=0:m-1
        s=s+v(j*2*n+i);
    end
    for j=0:m-1
        vHp(j*2*n+i)=s;
    end
end
vH=v-(1/m)*vHp;

end

