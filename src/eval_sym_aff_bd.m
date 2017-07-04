function [ min, max ] = eval_sym_aff_bd( coeff, const, s_min, s_max )
%EVAL_SYM_AFF_BD Summary of this function goes here
%   Evaluation bounds of symbolic expression: coeff * s + const
%   coeff: row vector with coefficients
%   const: constant
%   s: symbolic positive integer with bounds as input
%   s_min: column vector of minimum bounds
%   s_max: column vector of maximum bounds

% generate lower and upper bounds
n = length(s_max);
seq = 1:n;
temp = [coeff',seq'];    % compound matrix 
temp = sortrows(temp,1);	% sort according to lhs
nneg = sum(temp(:,1)<0);	% count negative entries in lhs
nzer = sum(temp(:,1)==0);	% count zero entries in lhs
npos = sum(temp(:,1)>0);	% count positive entries in lhs
% realign the matrix => entries in first column: pos, neg, zero:
temp = [temp(nneg+nzer+1:n,:);temp(1:nneg,:);temp(nneg+1:nneg+nzer,:)];

% calculate bounds
% MAX: max(positive terms) - min(negative terms) + const
max = temp(1:npos,1)'*s_max(temp(1:npos,2)) + temp(npos+1:npos+nneg,1)'*s_min(temp(npos+1:npos+nneg,2)) + const;
% MIN: min(positive terms) - max(negative terms) + const
min = temp(1:npos,1)'*s_min(temp(1:npos,2)) + temp(npos+1:npos+nneg,1)'*s_max(temp(npos+1:npos+nneg,2)) + const;

end

