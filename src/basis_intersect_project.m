function [ basis_bounds ] = basis_intersect_project( A, basis, b )
%BASIS_PROJECT 
%   project polytope (A<=b) to each basis vector 
%   return the lower and upper bounds on each vector

A_basis = A*basis;
disp('Number of basis vectors:');
disp(length(A_basis(1,:)));

basis_bounds = struct('lower', {}, 'upper', {});
for i=1:length(A_basis(1,:))
    str = sprintf('=== basis dim %d: ', i);
    disp(str);
    % Fourier-Motzkin Elimination modified to consider all final
    % inequalities
    % assume only b vector contains symbolic variables
    [lhs, rhs] = fourmotz(circshift(A_basis, -i, 2), b, 1);
    
    % generate lower and upper bounds
    m = length(lhs(:,1));
    temp = [lhs,rhs];			% compound matrix [lhs|rhs]
    temp = sortrows(temp,1);	% sort according to lhs
	nneg = sum(double(temp(:,1))<0);	% count negative entries in lhs
	nzer = sum(double(temp(:,1))==0);	% count zero entries in lhs
	npos = sum(double(temp(:,1))>0);	% count positive entries in lhs
    % realign the matrix => entries in first column: pos, neg, zero:
	temp = [temp(nneg+nzer+1:m,:);temp(1:nneg,:);temp(nneg+1:nneg+nzer,:)];

    % lower bounds (neg)
    lower = temp(npos+1:npos+nneg,2)./temp(npos+1:npos+nneg,1);
    disp('lower:');
    disp(lower);
    % upper bounds (pos)
    upper = temp(1:npos,2)./temp(1:npos,1);
    disp('upper:');
    disp(upper);

    basis_bounds(i).lower = lower; %ceil(lower);
    basis_bounds(i).upper = upper; %floor(upper);
end


end

