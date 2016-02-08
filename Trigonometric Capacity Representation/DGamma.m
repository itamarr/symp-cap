function [values] = DGamma(t,m,n,coeffs)
%% Gamma - we calculate the value of the derivative of Gamma at point t, coordinate wise. Gamma is represented by the coefficients in coeffs.
% For each coordinate of Gamma dot, we calculate its value at point t.
% Note that Gamma(t) = < coefficients, trigonometric par> 
% Here <,> is the standard scalar product.
    values = zeros(2*n, size(t,2));
    C = repmat(2*pi*[1:m],size(t,2),1).*cos([2*pi*t'*[1:m]]);
    S = -repmat(2*pi*[1:m],size(t,2),1).*sin([2*pi*t'*[1:m]]);
    trig = [cat(2, zeros(size(t')), S, C)];
    
    % It's probably possible to replace this for loop with one
    % multiplication of a repeated trig vector along the diagonal of a
    % matrix with the coefficients vector.
    for i=1:2*n

        values(i,:) = (trig * coeffs((i-1)*(2*m+1)+1: i*(2*m+1)))';
    end
end
