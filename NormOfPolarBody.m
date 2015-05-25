function [ n ] = NormOfPolarBody( P, y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [A, b] = PointsToHalfSpaces(P);
    A(~any(A,2),:) = [];
    b(~any(b,2),:) = [];
    options = optimset('Display','none'); %% disable logging for linprog
    maxX = linprog(-y,A,b,[],[],[],[],[],options);
    n = y * maxX;
end

