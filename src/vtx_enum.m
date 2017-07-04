function [ vtx ] = vtx_enum( A, b )
%VTX_ENUM Summary of this function goes here
%   Detailed explanation goes here

%% symbolic vertex enumeration : could be piece-wise parametric
% !!! assumption for simplicity: K is rectangular-like polyhedron 
% every pair of inequalities of parallel facets are ordered
% together in matrix A
disp('=================================================================');
disp('= Vertices Enumeration ');
disp('=================================================================');

n = length(A(1,:));

if(sum(abs(sum(A,2)) ~= ones(2*n,1) )|| sum(sum(A,1) ~= zeros(1,n)))
    error('Assumption is not satisfied: polytope of parametric iteration difference is not rectangular-like');
end

vtx = zeros(n, 2^n);
vtx = sym(vtx);
tmp = reshape(b,2,n);
tmp = diag([1,-1])*tmp;
% use binary representation to enumerate
for j = 0:2^n-1
%         str = sprintf('=== vertex %d ', j);
%         disp(str);
    v = zeros(n,1);
    v = sym(v);
    r = j;
    for k = n-1:-1:0
       v(n-k) = tmp(floor(r/(2^k))+1,n-k);
       r = mod(r,(2^k));
    end
    vtx(:,j+1) = v;
end
disp(vtx);


end

