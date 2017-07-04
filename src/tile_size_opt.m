function [ comp_time, cost_rec, solu_rec, mem_rec, ker_rec ] = tile_size_opt( ker, K_A, K_b, vtx, n, s, s_min, s_max, mem_tar, bnb_maxiter )
%TILE_SIZE_OPT Summary of this function goes here
%   tile size selection
%   using mixed integer geometric programming (YALMIP)
%   1) first MIGP with complete basis
%   2) remove redundant reuse basis according to previous tile size selection
%   3) iterate step 2 until tile size does not change

%% setup integer geometric programming
ns = length(s);
nk = length(ker);
t = sdpvar(ns,1); % t = s

% for mmm
% cond_unroll = t(2)<=31;
cond_init = [t>=s_min]; %

solu_rec = [];
mem_rec = [];
ker_rec = [];
cost_rec = [];
yalmip_rec = [];

disp('=================================================================');
disp('= First tile sizing with MIGP');
disp('=================================================================');

[solu_1, tm, tc, td, tSearch ]= solve_MIGP( ker, s, s_min, s_max, t, mem_tar, cond_init, s_max, bnb_maxiter );
solu_rec = [solu_rec;solu_1'];
mem_rec = [mem_rec;tm];
ker_rec = [ker_rec;ker];
cost_rec = [cost_rec;tc];
yalmip_rec = [yalmip_rec;td];

%% iterative opt
disp('=================================================================');
disp('= Iterative tile sizing');
disp('=================================================================');
solu_cur = solu_1;
solu_pre = s_max;
cond_redu = [];
while(any(solu_pre ~= solu_cur))  
    solu_pre = solu_cur;
    if(any(isnan(solu_pre)))
        error('Solver failed');
    end
    
    disp('#############################');
    disp('# Evaluate redundancy');
    disp('#############################');
    % evaluate redundancy
    redu = zeros(nk,1);
    for i = 1:nk % iterate over memory accesses
        str = sprintf('===== Kernal %d', i);
        disp(str);      
        redu(i) = 0;
        if(~isempty(ker(i).remove_basis)) 
            m = length(ker(i).remove_basis);
            seq = zeros(1, length(ker(i).reuse(1,:)));
            for j=1:m % iterate potential redundancy
                str = sprintf('=== Potential redundancy %d', j);
                disp(str);
                idx = isAlways(subs(ker(i).remove_basis(j).cond_ub, s, solu_pre) < 0); % be careful here !!!!!
%                 if (solu_pre(ker(i).remove_basis(j).num) == 0) % not helpful
%                     idx = ones(1,length(idx));
%                 end
                if( any(idx) )
                    % if previous solution makes redundancy exist           
                    disp('=== Found redundancy in reuse basis');
                    str = sprintf('* @ %d', ker(i).remove_basis(j).num);
                    disp(str);
                    disp('* cond:')
                    disp(ker(i).remove_basis(j).cond_ub(idx));
                    
                    % Handle redundancy
                    [A,b] = equationsToMatrix(ker(i).remove_basis(j).cond_ub(idx), s);
                    % collect redundancy condition
                    disp('* added constraint:')
                    str1 = sdisplay((double(A)*t)); % be careful here !!!!!
                    str2 = sdisplay(double(b-0.001));
                    disp(strcat(str1,' <= ',str2));
                    cond_redu = [cond_redu, (double(A)*t) <= double(b-0.001)]; % be careful here !!!!!

                    % record redundant reuse basis vectors
                    redu(i) = 1;
                    seq(ker(i).remove_basis(j).num) = 1;
                end
            end
            
            % remove corresponding reuse basis vector
            ker(i).reuse(:,any(seq,1)) = [];
            
        end
        
    end
    
    disp('#############################');
    disp('# Create new kernal models');
    disp('#############################');
    % Create new kernal model
    for i = 1:nk % iterate over memory accesses
        if(redu(i) == 1) 
            if(~isempty(ker(i).reuse))
                str = sprintf('===== Kernal %d', i);
                disp(str);  
                
                disp('===== Re-analyze redundancy of reuse basis ');
                disp(ker(i).reuse);
                [ ker(i).reuse, ker(i).remove_basis ] = proj_fme_redu( ker(i).reuse, K_A, K_b, s, s_min, s_max );

                disp('===== Create new full basis');
                [ ker(i).orth, ker(i).lat, ker(i).mem_basis] = create_full_basis( ker(i).reuse, vtx, n, s, s_max  );
            
            else
                error('ERROR: all reuse basis vectors are removed');
            end
        end
    end
    
    disp('#############################');
    disp('# new OPT problem');
    disp('#############################');
    % solve MIGP
    if(any(redu,1))
        disp('===== Solve new MIGP');
%         disp('* all added constraints:')
%         cond_redu
        [solu_cur, tm, tc, td] = solve_MIGP( ker, s, s_min, s_max, t, mem_tar, [cond_redu,cond_init], solu_pre, bnb_maxiter );
        solu_rec = [solu_rec;solu_cur'];
        mem_rec = [mem_rec;tm];
        ker_rec = [ker_rec;ker];
        cost_rec = [cost_rec;tc];
        yalmip_rec = [yalmip_rec;td];
    else
        disp('===== No redundancy found');
        solu_cur = solu_pre;
    end
    
end

% stop time record
comp_time = toc(tSearch);

disp('=================================================================');
disp('= Tile sizing finished');
disp('=================================================================');

disp('=================================================================');
disp('= Optimization record:');
disp('=================================================================');

disp('=== computation time:  ');
disp(comp_time);

disp('=== solutions');
disp('* max:');
disp(s_max');
disp('* records:');
disp(solu_rec);

disp('=== yalmip info');
for i = 1:length(yalmip_rec)
    str = sprintf('* Opt %d: ', i);
    str = strcat(str, yalmiperror(yalmip_rec(i)));
    disp(str);
end
disp(' ');

disp('=== tile size');
disp(prod(solu_rec+1,2))

disp('=== mem size');
str = sprintf('* target : %d', mem_tar);
disp(str);
disp('* records:');
disp(mem_rec);

disp('=== tile access cost');
disp(cost_rec);

% return minimum
[cost_rec, idx] = min(cost_rec);
solu_rec = solu_rec(idx,:);
mem_rec = mem_rec(idx);
ker_rec = ker_rec(idx,:);
if(idx>1)
    disp('!!!!!!! redundancy help to find better solution');
end
    
% clear yalmip variables
yalmip('clear');

end

