function [ T ] = rowcol_add(  q, dim, idx1, idx2 , n )
%ROWCOL_ADD Summary of this function goes here
%   Detailed explanation goes here
%   dim == 1: add row(idx2) with -q times row(idx1)
%   A = T*A
%   dim == 2: add col(idx2) with -q times col(idx1):
%   A = A*T
%   return the tansformation matrix T

% n is the length of row/col
T = eye(n);

if(dim == 1)
    % row 
    T(idx2,idx1) = -q;
    T(idx2,idx2) = 1;    
elseif(dim == 2)
    % column
    T(idx1,idx2) = -q;
    T(idx2,idx2) = 1; 
end

end

