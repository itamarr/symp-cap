function [ grad ] = GradOfPolarNorm( P, y )
    [M,I] = max(P*y); %vertices of P correspond to normal vectors of P^o,
                      %Find which has smallest angle w.r.t y.
    grad = P(I(1),:);
  
end