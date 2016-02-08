function [y] = dF(coeff,P,m,n)

mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)];
y=zeros(size(coeff));
for k = 1:m
    for j = 1:n
        y((j-1)*(2*m+1)+1+k) = integral(@(x) integrandA(mJ2n, coeff, m, n, k, j, x, P), 0, 1); %d/da
        y((j-1)*(2*m+1)+1+m+k) = integral(@(x) integrandB(mJ2n, coeff, m, n, k, j, x, P), 0, 1); %d/db
        %d/dconst = 0  gammadot has 0 there
    end
end


end

function [result] = integrandA(mJ2n, coeff, m, n, k, j, x, P)
    gradG = (dG(mJ2n * DGamma(x, m, n, coeff), P, n));
    v = zeros(2*n,1);
    v(j) = 1;   %% construct stdrd basis vect, apply J to it, repeat and take product with each column of grad G
    v = mJ2n*v;
    repV=repmat(v,1,size(x,2));
    result = (dot(gradG,repV)).*(-2*pi*k*sin(2*pi*k*x));
end

function [result] = integrandB(mJ2n, coeff, m, n, k, j, x, P)
    gradG = (dG(mJ2n * DGamma(x, m, n, coeff), P, n));
    v = zeros(2*n,1);
    v(j) = 1;
    v = mJ2n*v;
    repV=repmat(v,1,size(x,2));
    result = (dot(gradG,repV)).*(2*pi*k*cos(2*pi*k*x));
end