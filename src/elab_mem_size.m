function [ mem_all, mem_all_real, mem_b_real, mem_all_sym, cost, cost_global, cost_local, cost_sym, cost_min ] = elab_mem_size(  ker, s, s_max, t )
%elab_mem_size Summary of this function goes here
%   Elaborate the expression of memory size and communication cost
% 

ns = length(s);
nk = length(ker);

disp('=================================================================');
disp('= Elaborate the expression of memory size and communication cost');
disp('=================================================================');
mem_b = sdpvar(nk+1,1);
mem_b_real = sdpvar(nk+1,1);
mem_b_max = ones(nk,1);
mem_b_sym = sym(ones(nk,1));

% l = [];
% cond_lat = [];
for i = 1:nk % iterate over memory accesses
    str = sprintf('=== ker(%d) ', i);
    disp(str);
    nl = length(ker(i).lat);
    
    %% foreach diminsion of the full basis matrix
    %disp('* basis memory:')
    mem_b(i,1) = 1;
    mem_b_real(i,1) = 1;
    for j=1:length(ker(i).mem_basis(:,1))
        % determine expression
        if(j<=(ns-nl))
            %% reuse basis vectors
            if (ker(i).mem_basis(j,2) == 1)
                continue;
            else
                sw = (double(subs(ker(i).mem_basis(j,1),s,s_max)) <= double(ker(i).mem_basis(j,2)) ); 
                if(sw)
                    % expr from rows in G
                    [A,b] = equationsToMatrix(ker(i).mem_basis(j,1), s);             
                else 
                    % expr from diagonal elements in S
                    [A,b] = equationsToMatrix(ker(i).mem_basis(j,2), s);
                end
            end
        else
            %% orth basis vectors
            % apply lattice
            str = sprintf('!!! Apply lattice factor on dim(%d) ', j);
            disp(str);
            
            sw = 0;
            nu = length(ker(i).lat(j-ns+nl).upper(:,1)); % number of upper bounds
            for k=1:nu
                % Compare: diff(G(j,:)*x) and s(j)*(lat_upper(j) + epsilon)
                % epsilon represents the round up error/distance in [0~1] 
                upper = ker(i).lat(j-ns+nl).upper(k,1:ns)*s + ker(i).lat(j-ns+nl).upper(k,ns+1);
                frac = any(ker(i).lat(j-ns+nl).upper(k,1:ns) < 1, 1) - any(ker(i).lat(j-ns+nl).upper(k,1:ns) == 0, 1);
                frac = any(frac>0);
                % check delta = diff(G(j,:)*x) - s(j)*lat_upper(j)
                delta = ker(i).mem_basis(j,1) - ker(i).mem_basis(j,2) * upper;
                [A,b] = equationsToMatrix(delta, s);
                if (all([A,-b]>=0, 2))
                    % when delta - s(j) * epsilon < 0
                    %   => epsilon > max( delta/s(j) )
                    %disp('!!! always positive delta')
                    eps_min = (double(A)*s_max - double(b)) / double(ker(i).mem_basis(j,2));
                    if (eps_min < 0.1 && frac || eps_min == 1 && ~frac)
                        % diff(G(j,:)*x) is expected to be smaller mostly
                        % or equal
                        %disp('!!! delta is relatively much smaller than epsilon')
                        sw = 1;
                        break;
                    end
                elseif (all(A<=0, 2))
                    if (double(A)*s_max - double(b) < 0)
                        %disp('!!! always negative delta')
                    	% diff(G(j,:)*x) is sure to be smaller
                        sw = 1;
                        break;
                    end
                end    
            end
            
            % determine expression
            if(sw)
                % expr from rows in G
                disp('!!! use expr from G')
                [A,b] = equationsToMatrix(ker(i).mem_basis(j,1), s);  
            else
                % expr from diagonal elements in S*L
                disp('!!! use expr from diag(S*L)')
                upper = ker(i).lat(j-ns+nl).upper(:,1:ns)*s + ker(i).lat(j-ns+nl).upper(:,ns+1);
                [A,b] = equationsToMatrix(ker(i).mem_basis(j,2)*upper, s);
            end
            
        end    

        %% determine expression
        expr_max = double(A)*s_max-double(b);
        expr_real = double(A)*t-double(b);
        expr = expr_real;
        
        expr_sym = A*s-b;
        
        % apply floor() for lattice-applied dimensions
        if (j>(ns-nl) && sw == 0)
            [~,den] = numden(sym(A));
            den = max(double(den),[],2);
            idx = any(den>1,2);
            expr_real(idx) = floor(expr_real(idx));
            expr_max(idx) = floor(expr_max(idx));
            expr_real = expr_real + 1;
            expr_max = expr_max + 1;
            expr = expr + 1;
            
            expr_sym(idx) = floor(expr_sym(idx));
            expr_sym = expr_sym + 1;
        end
        
        %% calculate memory size
        %disp(sdisplay2(expr)) % NOTE: affect speed in RandEnum test 
        mem_b(i,1) = mem_b(i,1) * max(expr); % posynomial
        mem_b_real(i,1) = mem_b_real(i,1) * max(expr_real); % posynomial and include floor()
        mem_b_max(i,1) = mem_b_max(i,1) * max(expr_max);
        
        % assume no need to max
        if(size(A,1)>1)
           error('Need to apply max()');
        end
        mem_b_sym(i,1) = mem_b_sym(i,1) * expr_sym; 
    end

end

% for single kernal case
mem_b(nk+1,1) = 0;
mem_b_real(nk+1,1) = 0;

%% total memory requirements
mem_all = sum(mem_b,1); %norm 1
mem_all_real = sum(mem_b_real,1); 
mem_all_sym = sum(mem_b_sym, 1);

disp(' ')
disp('=== Max memory usage:')
disp(sum(mem_b_max,1));
disp('== Memory usage breakdown:')
disp(mem_b_max');

% find significantly small arrays
[~, idx] = max(mem_b_max);
mem_scale = mem_b_max(idx) ./ mem_b_max;
mem_small = mem_scale>10;

%% total communication cost
cost = 0;
cost_global = 0;
cost_sym = 0;
cost_min = 0;
for i = 1:nk 
    % real cost
    cost = cost + ker(i).coeff * mem_b_real(i,1) * prod((s_max(1:ker(i).inner_pos)+1)./(t(1:ker(i).inner_pos)+1));
    
    % cost for GP solver
    % feasible with geometric programming 
    % when variables are relaxed as continuous
    cost_global = cost_global + ker(i).coeff * mem_b(i,1) * prod(s_max(1:ker(i).inner_pos)+1) / prod(t(1:ker(i).inner_pos));
      
    % symbolic cost
    cost_sym = cost_sym + ker(i).coeff * mem_b_sym(i,1) * prod((s_max(1:ker(i).inner_pos)+1)./(s(1:ker(i).inner_pos)+1));
    
    % minimum cost
    mini_coeff = 1;
    if ker(i).rw == 3
        mini_coeff = 2;
    end
    cost_min = cost_min + mini_coeff * mem_b_max(i,1);
    
end

%% cost for local
% feasible with BNB/nonlinear-fmincon
% cost_local = 0;
% for i = 1:nk            
    % 1) fast but not close to optimum
%     cost_local = cost_local - ker(i).coeff * prod(t(1:ker(i).inner_pos)+1);
    % 2) slow but quite close to optimum
%     cost_local = cost_local * ker(i).coeff * prod(t(1:ker(i).inner_pos)+1);        
% end

coeff = zeros(ns,nk);
for i = 1:nk
    coeff(1:ker(i).inner_pos,i) = 1;
end
coeff(:,mem_small) = [];
coeff = unique(coeff', 'rows');
coeff = sum(coeff', 2);

% cost_local = -prod(t+1); %OLD: transfer per tile
cost_local = -prod(t.^coeff+1); % find different results on macOS and Windows
% cost_local = -prod(coeff.*(t.^coeff)+1); % most balanced for cnnlayer1



end



