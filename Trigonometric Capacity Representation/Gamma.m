function [values] = Gamma(t,m,n,coeffs)
%% Gamma - we calculate the value of Gamma at point t, coordinate wise. Gamma is represented by the coefficients in coeffs.
% For each coordinate of Gamma, we calculate its value at point t.
% Note that Gamma(t) = < coefficients, trigonometric par> 
% Here <,> is the standard scalar product.
    values = zeros(1,2*n);
    C = cos([2*pi*[1:m]*t]);
    S = sin([2*pi*[1:m]*t]);
    trig = [1 cat(2, C, S)];
    
    % It's probably possible to replace this for loop with one
    % multiplication of a repeated trig vector along the diagonal of a
    % matrix with the coefficients vector.
    for i=1:2*n

        values(i) = trig * coeffs((i-1)*(2*m+1)+1: i*(2*m+1))
    end

end
