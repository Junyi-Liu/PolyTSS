%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Memory mapping for Data Reuse 
% with a Parametric Iteration difference Space
% Application: tile size selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: symbolic/parametric Polyhedra in MPT is not supported

% parametric_reuse_mapping
function [ ker ] = PRM( A, b, vtx, F, rw, s, s_min, s_max)

n = length(A(1,:));

%% Define result object: ker
% reuse:    basis vectors of reuse kernal
% orth:     basis vectors of orthogonal kernal space to reuse kernal
% cond:     all(cond < 0), conditions of parameters for current resue basis
%           without redundancy
% lat:      diagonal elements of the rectangular lattice on orth space
%           (only intersect iteration difference polytope proejection 
%            at 0 point)
%           Full expression: floor(min(l(x).upper) + 1); 

% NOTE: 1) the first kernal has the orignal reuse basis without removing
%           redundancy
%       2) the final memory requirement = det(basis) * det(diag(l))

ker = struct('reuse', {}, 'orth', {}, 'remove_basis', {}, ...
             'G', {}, 'diagS', {}, 'lat', {}, 'mem_basis', {}, ...
             'rw', {}, 'coeff', {}, 'inner_pos', {}, 'all_pos', {} );
ker(1).rw = rw;

% position of innermost iterators
ker(1).all_pos = any(F,1);
[~, reverse_pos] = max(fliplr(ker(1).all_pos));
ker(1).inner_pos = (n+1) - reverse_pos; 

% coefficient
ker(1).coeff = 1;
if ker(1).rw == 2 % for write access
    % related innermost loop position is larger than the number of related iterators
    if ker(1).inner_pos > sum(ker(1).all_pos) 
        ker(1).coeff = 2;
    end
end
if ker(1).rw == 3 % for read-and-write access
    ker(1).coeff = 2;
end

%% parametric iteration difference set: 0-symmetric polytope K
%  (Minkowski addition of iteration space D and -D)
%  A*v <= b, 0-symmetric polytope 

%% array access reference: Fv+c, data reuse charaterization 
% find the basis of the nullspace of F

disp('=================================================================');
disp('= Find the kernal basis of the array reference ');
disp('=================================================================');

% % FCCM07 example: A[50*v0 + v1 + v2]
% F = [50 ,1 ,1];
% q = 1;
% n = 3;

% FCCM07 full loop: A[4*m + 50*v0 + v1 + v2]
%F = [4, 50 ,1 ,1];
q = length(F(:,1));

% Hermite Normal Form
% F*V' = H'
[V, H] = hermiteForm(F');
% replace by self-transpose 
H = H';
V = V';

% Use unimodular matrix V to form the basis
ker_reuse = V * [zeros(q,1), zeros(q,n-q); zeros(n-q,1), eye(n-q)];
ker_reuse(:, ~any(ker_reuse,1)) = []; % eliminate zero column vectors

% FCCM07 example 
%ker_reuse = [0; 1; -1];

disp(ker_reuse);

%% Remove redundant basis vectors by projection
disp('=================================================================');
disp('= Projection of polytope K intersecting the kernal subspace');
disp('=================================================================');

[ ker_reuse, remove_basis ] = proj_fme_redu( ker_reuse, A, b, s, s_min, s_max );

% complete reuse basis
ker(1).reuse = ker_reuse; 

% add conditions
if(~isempty(remove_basis))
    disp('###### found conditions of removing redundant basis vectors');
    ker(1).remove_basis = remove_basis;
end

% when 2nd basis vector can be removed
% ker_reuse(:,2) = [];
% disp(ker_reuse);

disp('*****************************************************************');
disp('* Reuse basis with all possible redundancy')
disp('*****************************************************************');


if(isempty(ker_reuse))
    warning('no reuse basis found');
    ker(1).orth = [];
    ker(1).lat = [];
    ker(1).mem_basis = sym(ones(n,2));    
    for k=1:length(F(:,1))
        ker(1).mem_basis(k,1) = F(k,:)*s+1;
        ker(1).mem_basis(k,2) = inf;        
    end
else
    [ ker(1).orth, ker(1).lat, ker(1).mem_basis, ker(1).G, ker(1).diagS] = create_full_basis( ker(1).reuse, vtx, n, s, s_max  );
end

disp('=================================================================');
disp('= FINISH');
disp('=================================================================');






