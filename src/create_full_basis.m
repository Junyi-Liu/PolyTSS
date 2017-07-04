function [ orth, proj_orth, mem_basis, G, diagS] = create_full_basis( reuse, vtx, n, s, s_max  )
%CREATE_FULL_BASIS Summary of this function goes here
%   find orth basis according to input resue basis
%   project vertices onto orth space to determine rectangular lattice
%   generate full basis and its memory size


%% Find Orthogonal subspace to reuse space

disp('=================================================================');
disp('= Orthogonal basis to the reuse kernal');
disp('=================================================================');

% Hermite Normal Form
m = length(reuse(1,:));
Q = reuse';
[V, H] = hermiteForm(Q');
H = H';
V = V';

% Use unimodular matrix V to form the basis
orth = V * [zeros(m,1), zeros(m,n-m); zeros(n-m,1), eye(n-m)];
orth(:, ~any(orth,1)) = []; % eliminate zero column vectors

disp(orth);

%% K projection to the orthogonal subspace
disp('=================================================================');
disp('= Project iteration difference set to ker_orth space');
disp('=================================================================');

% project vertices to the subspace of orth basis
disp('====== Project vertices to calculate coordinates in orth space');
T = eye(n-m)/(orth'*orth) * orth';
vtx_orth = T*vtx;

% vtx_polyhedron_plot2D( vtx_orth, s, n, 100 );

%% Find dense lattice in the subspace orthogonal to reuse basis

disp('=================================================================');
disp('= Find dense lattice in the subspace orthogonal to reuse basis');
disp('=================================================================');

% Intuition: rectangular lattice
% create diagonal elements of dense lattice
% The projection of polytope K on the orthogonal space is still
%   0-symmetric
% Full expression: floor( max(proj_orth(x).upper) + 1 )
% Note: use floor with +1, since max(proj_orth(x).upper) could be integer
disp('=== rectangular lattice');

% take the maximum positive value on each axis 
% Note: 1) we take all positive values for symbolic execution
%       2) the projected polytope is also 0-symmetric
disp('====== Check projection bounds on each space dimensions');
proj_orth = struct('upper', {});
if(~isempty(s))
    for j = 1:n-m
        str = sprintf('=== Orth space dim: %d ', j);
        disp(str);
        if(any(has(vtx_orth(j,:), s)))
            % T*s = c
            [T,c] = equationsToMatrix(vtx_orth(j,:),s);
            ns = length(s);
            % remove identical rows
            T = unique([T,-c], 'rows');
            % remove zero points
            T(any(all(double(T) == 0, 2)), :) = [];
            %T
            % filter points
            if(any(all(double(T) >= 0, 2)))
                % keep always positive value
                disp('** Found always positive value');
                T = T(all(double(T) >= 0, 2),:);
                seq = (double(T(:,1:ns)*s_max + T(:,ns+1))<1);
                if(all(seq,1))
                    %proj_orth(j).upper = 0;
                    proj_orth(j).upper = zeros(1,ns+1);
                else
                    %proj_orth(j).upper = T(~seq,1:ns)*s + T(~seq,ns+1);
                    proj_orth(j).upper = double(T(~seq,:));
                end
            elseif(isempty(T))
                % just zero point
                proj_orth(j).upper = zeros(1,ns+1);
            else
                % remove always negative value
                T(all(double(T) <= 0, 2),:) = [];
                %proj_orth(j).upper = T(:,1:ns)*s + T(:,ns+1);
                proj_orth(j).upper = double(T);
            end
            disp('** coefficient matrix of symbolic variables:')
            disp(proj_orth(j).upper);
            disp('** positive upper bounds:');
            disp(proj_orth(j).upper(:,1:ns)*s + proj_orth(j).upper(:,ns+1)); 
        else
            % keep all potentially positive values
%             proj_orth(j).upper = vtx_orth(j, ~any(vtx_orth(j,:) < 0, 1));
%             [row,col] = size(proj_orth(j).upper);
%             proj_orth(j).upper = reshape(proj_orth(j).upper, [col,row]);
            proj_orth(j).upper = max(vtx_orth(j, ~any(vtx_orth(j,:) < 0, 1)));
            disp('** positive upper bound:');
            disp(proj_orth(j).upper); 
        end
    end
else
    disp('=== No symbolic variables');
    %mem_orth = 1;
    for j = 1:n-m
        str = sprintf('=== Orth space dim: %d ', j);
        disp(str);
        proj_orth(j).upper = max(vtx_orth(j, ~any(vtx_orth(j,:) < 0, 1)));
        disp(proj_orth(j).upper);
        %mem_orth = mem_orth * floor(proj_orth(j).upper+1);
    end      
end


%% Full basis and its Smith/Diagonalized Form

disp('=================================================================');
disp('= Full basis and its Smith/Diagonalized Form');
disp('=================================================================');

ker_full = [reuse,orth];
disp('Full basis without scaling on orthogonal space: ');
disp(ker_full);

% smith form: S = G*A*U
%     [G, U, S] = smithForm(ker_full);

% diagnoal form
[G, U, S] = diagonalizer(ker_full, n, n);
G
U
S
% idx

% verify (L^-1 * U * L) is still a unimodular matrix
for i = 1 : n-m
    tmp = U(i+m,:);
    tmp(i+m) = 0;
    if(any(tmp,2))
        idx = 1:n;
        idx = idx(any(tmp,1));
        str = sprintf('%d, ', idx);
        str = sprintf('!!! please check l(%d) can be dividied by l(%d)', str, i);
        disp(str);
        warning('(L^-1 * U * L) may not be unimodular' )
        error('Need verification on divisability between lattice diagonal elements');
    end
end

% process G and S for memory calculation
diagS = double( diag(S) );
neg_mask = ones(size(diagS));
neg_mask(diagS<0) = -1;
diagS = abs(diagS);
if(~isempty(s))
    mem_basis = sym(ones(n,2));
    
    % full basis
    for k=1:n
         % negate all coefficients due to negative modular value
        G(k,:) = neg_mask(k)*G(k,:); 
        % make the negative cofficient positive for diff =
        % max(G*x)-min(G*x), where x in IterDom, s in IterDiffDom (could be negative)
        % a straight way to calculate the maximum interval on current
        % dimension 
        %G(k,G(k,:)<0) = G(k,G(k,:)<0)*(-1);
        tmp = abs(G(k,:));
        % assign the interval expression and modular value
        mem_basis(k,1) = tmp*s+1;
        mem_basis(k,2) = diagS(k);        
    end
    
else
    mem_basis = 1;
    for k=1:n
        if(S(k,k) > (G(k,:)*s_max + 1))
            mem_basis = mem_basis * (neg_mask(k)*G(k,:)*s_max+1);
        else
            mem_basis = mem_basis * abs(S(k,k));
        end
    end
end
    
    

end


