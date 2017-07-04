function [ T ] = rowcol_itrchg( idx1, idx2 , n)
%ROWCOL_ITRCHG Summary of this function goes here
%   Detailed explanation goes here
%   interchange row/col(idx1) and row/col(idx2) of a matrix A
%   return the tansformation matrix T
%   row interchange: A = T*A 
%   col interchange: A = A*T 

% n is the length of row/col
T = eye(n);

% row/col idx1
T(idx1,idx1) = 0;
T(idx1,idx2) = 1;

% row/col idx2
T(idx2,idx2) = 0;
T(idx2,idx1) = 1;

end

