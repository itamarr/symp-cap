function [y] = F(coeff,P,m,n)

mJ2n = [zeros(n),-eye(n);eye(n),zeros(n)];



y = integral(@(x) G(mJ2n*DGamma(x, m, n, coeff), P, n), 0, 1);

end