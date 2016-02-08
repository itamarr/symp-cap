%n =1
%m=7
%coeffs = rand(18,1);
results = zeros(2, 200);
for t=1:200
    results(:, t) = Gamma(t/200, m,n, gamma);
end
scatter(results(1,:), results(2,:))
%{
m = 2;
n = 2;
sz = 2*n*(2*m+1);
B = zeros(sz,sz);
coeffs = [0 1 2 1 1 0 1 2 2 2 0 1 2 3 3 0 1 2 4 4];
for k = 1:m
    for j = 1:n
        if (((n+j-1)*(2*m+1)+1+k > sz) || ((j-1)*(2*m+1)+1+m+k > sz) || ((j-1)*(2*m+1)+1+k > sz) || ((n+j-1)*(2*m+1)+1+m+k > sz))
            print "wtf"
        end
        B((n+j-1)*(2*m+1)+1+k, (j-1)*(2*m+1)+1+m+k) = k;
        B((j-1)*(2*m+1)+1+k, (n+j-1)*(2*m+1)+1+m+k) = -k;
    end
end

%}