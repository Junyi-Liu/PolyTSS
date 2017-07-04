function [ comp_time, cost_rec, solu_rec, mem_rec, coverage_rec ] = tile_size_enum_rand( ker, s, s_min, s_max, mem_tar, t_enum )
%TILE_SIZE_ENUM Summary of this function goes here
%   tile size selection by enumeration
%   

%% setup integer geometric programming
ns = length(s);
nk = length(ker);

%% mem size
t = sdpvar(ns,1); 
[ ~, ~, ~, mem_all_sym, ~, ~, ~, cost_sym, ~ ] = elab_mem_size(  ker, s, s_max, t );

% start time record
tic;

% % symbolic access cost
% % cost = comm_all_real * prod((s_max+1)./(t+1));
% cost = comm_all_sym * prod((s_max+1)./(s+1));

%% enumeration
disp('=================================================================');
disp('= Tile size random enumeration');
disp('=================================================================');

% control initialization
mem_cur = mem_tar;
s_cur = zeros(ns,1);

% random number generator setting
% state_rng = rng(0); 
state_rng = rng('shuffle');

% record initialization
tm = inf;
tc = inf;
ts = zeros(ns,1);
cnt = 0;

% random search
nsp = prod(s_max - s_min +1);
str = sprintf('===== number of all candidates: %d ', nsp);
disp(str);

% search loop 
comp_time = 0;
while(comp_time <= t_enum) %mem_cur <= mem_tar
    
    % randomly update tile sizes
    % Uniformly distributed pseudorandom integers
    for i=1:ns
        s_cur(i) = randi([s_min(i), s_max(i)]);
    end
    
    % show progress
%     str = sprintf('=== Progress: ( %d / %d ) seconds', comp_time, t_enum);
%     disp(str);
    cnt = cnt+1;
%     disp(s_cur');
    
    % check candidate
%     assign(t, s_cur);
%     mem_cur = value(mem_all_real);
    mem_cur = subs(mem_all_sym, s, s_cur);
    if(mem_cur <= mem_tar) 
%         cost_cur = value(cost);
        cost_cur = subs(cost_sym, s, s_cur);
        if (cost_cur < tc) 
%             disp('# found better solution')
            tc = double(cost_cur);
            ts = s_cur;
            tm = double(mem_cur);
%             disp('* solution');
%             disp(ts');
%             disp('* mem');
%             disp(tm);
%             disp('* cost');
%             disp(double(tc));
        else
%             disp('# skip solution: larger cost')
        end
    end
    
    % current elapsed time 
    comp_time = toc;
    
end

% stop time record
comp_time = toc;

solu_rec = ts';
mem_rec = tm;
cost_rec = tc;
coverage_rec = cnt/nsp;

disp('=================================================================');
disp('= Tile size enumeration finished');
disp('=================================================================');

% diary on

disp('=================================================================');
disp('= Enumeration result:');
disp('=================================================================');

disp('=== computation time:  ');
disp(comp_time);

% number of random searchs 
disp('=== number of random search: ');
disp(cnt);

disp('=== solutions');
disp('* max:');
disp(s_max');
disp('* records:');
disp(solu_rec);

disp('=== mem size');
str = sprintf('* target : %d', mem_tar);
disp(str);
disp('* records:');
disp(mem_rec);

disp('=== tile access cost');
disp(cost_rec);

% diary off

% clear yalmip variables
yalmip('clear');

end

