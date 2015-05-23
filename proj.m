function [zH] = proj(z,m,n)
%PROJ Projection onto gi=0-plane
%   Projection of a vector onto the H-plane
%   described in section 2.3 of the paper, equation (2.17)

zHp=zeros(2*m*n,1); %the sum in the brackets

for i=1:2*n %go through the index i of the sum
    
    s=0; %will be the scalar product
    %calculating the scalar product
    for j=0:m-1
        s=s+z(j*2*n+i);
    end
    
    %Multiplying the scalar product with a(i)
    % which is equivalent to filling up certain components of zHp
    for j=0:m-1
        zHp(j*2*n+i)=s;
    end
    
end

%computing z as in the equation (2.17)
zH=z-(1/m)*zHp;

end

