function [ ker_reuse, remove_basis ] = proj_fme_redu( ker_reuse, A, b, s, s_min, s_max )
%PROJ_FME_REDU Summary of this function goes here
%   Detailed explanation goes here
%   ker_reuse: reuse basis of access Kernal
%   remove_basis: record conditions of potential redundant basis vectors

o = 1;
while(o>0)
    o = 0;  
    disp('###### FME Projection');
    m = length(ker_reuse(1,:));
    proj_reuse = basis_intersect_project(A, ker_reuse, b);

    % A basis vector is redundant when its projection :
    % max(lower bound) > -1 and min(upper bound) < 1
    % => remove condition: max(cond_lb) > -1 && min(cond_ub) < 1
    disp('###### Find Redundancy');
    lb = struct('keep', {}, 'remove', {});
    ub = struct('keep', {}, 'remove', {});
    seq = zeros(1,m);
    for i = 1:m
        str = sprintf('=== basis vector(%d) ', i);
        disp(str);
        
        % lower bound vector
        %     lb_keep = isAlways(proj_reuse(i).lower <= -1);
        %     lb_remove = isAlways(proj_reuse(i).lower > -1);
        if(isempty(s))
            lb(i).keep = (proj_reuse(i).lower <= -1);
            lb(i).remove = (proj_reuse(i).lower > -1);
        else
            [TA,Tb] = equationsToMatrix(proj_reuse(i).lower,s);
            lb_num = length(TA(:,1));
            lb(i).keep = zeros(lb_num, 1);
            lb(i).remove = zeros(lb_num, 1);
            for j = 1:lb_num
                [lb_min,lb_max] = eval_sym_aff_bd(double(TA(j,:)), -double(Tb(j)), s_min, s_max);
                lb(i).keep(j) = (lb_max<=-1);
                lb(i).remove(j) = (lb_min > -1);
            end
        end
        % upper bound vector
        %     ub_keep = isAlways(proj_reuse(i).upper >= 1);
        %     ub_remove = isAlways(proj_reuse(i).upper < 1);
        if(isempty(s))
            ub(i).keep = (proj_reuse(i).upper >= 1);
            ub(i).remove = (proj_reuse(i).upper < 1);
        else
            [TA,Tb] = equationsToMatrix(proj_reuse(i).upper,s);
            ub_num = length(TA(:,1));
            ub(i).keep = zeros(ub_num, 1);
            ub(i).remove = zeros(ub_num, 1);
            for j = 1:ub_num
                [ub_min,ub_max] = eval_sym_aff_bd(double(TA(j,:)), -double(Tb(j)), s_min, s_max);
                ub(i).keep(j) = (ub_min >= 1);
                ub(i).remove(j) = (ub_max < 1);
            end
        end
        
        % Find redundant basis vectors
        if(any(lb(i).remove) && any(ub(i).remove))
            % basis can be removed certainly
            disp('removed');
            seq(1,i) = 1;
            % record basis change
            o =o+1;
%             ker_reuse(:,i) = [];
%             % basis index offset
        end
    end
    
    % remove determined redundant basis vectors
    ker_reuse(:,any(seq,1)) = [];
    disp('###### Current Reuse Kernal');
    disp(ker_reuse);
    
    % record redundancy when no determined basis change
    if(o==0)
        disp('###### Record Potential Redundancy');
        k = 1;
        remove_basis = struct('num', {}, 'cond_ub', {});
        for i = 1:m
            % create conditions
            str = sprintf('=== basis vector(%d) ', i);
            disp(str);
            if(all(lb(i).keep,1) && all(ub(i).keep,1))
                % basis must be kept
                disp('kept');
                continue;
            else
                % find conditions for potential redundancy
                expr_list = proj_reuse(i).upper(~any(ub(i).keep,2));       
                % eliminate denominator
                [~,den] = numden(expr_list);

                r = zeros(1, length(expr_list(:,1)));
                for j = 1:length(expr_list(:,1))               
                    % expr < 1 => expr*lcm - lcm < 0
                    tmp = lcm(den(j,:));
                    expr_list(j,:) = tmp.*expr_list(j,:) - tmp;
                    % record constraints for " expr < 1 "
                    if(tmp==1)  
                        r(j) = 1;
                    end
                end                

                % drop recorded constraints  
                expr_list(any(r,1),:) = [];
                
                % record symbolic expressions
                if(~isempty(expr_list))                    
                    remove_basis(k).cond_ub = expr_list;
                    disp('potentially redundant');
                    remove_basis(k).num = i;
                    disp('condition:');
                    disp(remove_basis(k).cond_ub);
                    k = k+1;
                end
            end
        end
    end
    
end
str = sprintf('###### number of possible redundant basis vectors : %d ', k-1);
disp(str);

end

