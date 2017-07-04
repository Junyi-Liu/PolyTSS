function [ ts, tm, tc, td, tSearch ] = solve_MIGP( ker, s, s_min, s_max, t, mem_tar, cond_init, t_init, bnb_maxiter )
%SOLVE_MIGP Summary of this function goes here
%   solve specific MIGP problem
%   variable t = s, diff tile size

ns = length(s);
    
%% mem size
[ mem_all, mem_all_real, mem_b_real, ~, cost, cost_global, cost_local, ~, cost_min ] = elab_mem_size(  ker, s, s_max, t );

% start time record
tSearch = tic;

%% solver settings
ops = sdpsettings;

% control display
ops.verbose = 0;
ops.warning = 0;
ops.showprogress = 0;
% ops.debug = 1;

% control algorithm
% ops.fmincon.Algorithm = 'sqp';
% ops.fmincon.Algorithm = 'active-set';

% control start point
% ops.usex0 = 1;
% assign(t,t_init);

%% solve optimization problem
disp('=================================================================');
disp('= GP solver run');
disp('=================================================================');

% geometric programming for global minimum
% NOTE: 1) variables are relaxed to be continuous 
%       2) local minimum (Matlab default)
ops.relax = 0;
ops.solver = 'fmincon-geometric';
% ops.fmincon.Algorithm = 'sqp'; 
% ops.usex0 = 1;
% assign(t, s_max);

obj = cost_global;
cond = [t>=0, t<=s_max, mem_all_real<=mem_tar, cond_init]; %cond_lat, l>=1, 
% cond = [t>=0, t<=s_max, mem_all<=mem_tar, cond_init]; 

diagnostics = optimize(cond, obj, ops); %
disp(diagnostics);

if diagnostics.problem == 3
%     str = yalmiperror(diagnostics.problem);
%     warning(str);
elseif diagnostics.problem == 1
     % infeasible
    str = yalmiperror(diagnostics.problem);
    warning(str);
    assign(t, zeros(ns,1));
elseif diagnostics.problem ~= 0
    str = yalmiperror(diagnostics.problem);
    error(str);
else
%     str = yalmiperror(diagnostics.problem);
%     disp(str);
end

% evaluate solution
t_GP = round(value(t)); 
% sw = any(floor(value(t)) == 0); % if some variables tends to be zeros.
sw = any(value(t) <= 1.5); % if some variables tends to be zeros.
t_zero = any(value(t) < 1, 2);

% simple approximation by rounding
disp('=== Continuous tile size selection :');
disp(value(t'));
disp('=== Rounded tile size selection :');
disp(t_GP');

assign(t, t_GP);
cost_GP = value(cost);
mem_GP = value(mem_all_real);
disp('=== Memory:');
disp(value(mem_all_real))
disp('=== tile access cost:');
disp(value(cost))

disp('=================================================================');
disp('= MIGP/non-linear BNB solver run');
disp('=================================================================');
% BNB with non-linear solver
% NOTE: local minimum

% config
if (~isempty(bnb_maxiter))
    ops.bnb.maxiter = bnb_maxiter;
end
ops.relax = 0;
ops.solver = 'bnb';
% ops.bnb.solver = 'fmincon-geometric';
% ops.fmincon.Algorithm = 'trust-region-reflective'; 

ops.usex0 = 1;
assign(t, t_GP);

% sw = 1;
if (sw)
    disp('===== use local obj')
%     obj = cost_local;
%     ops.fmincon.Algorithm = 'trust-region-reflective';

    obj = cost_local;
    
else
    disp('===== use global obj')
    obj = cost_global;   
end

cond = [integer(t), t>=0, t<=s_max, mem_all_real<=mem_tar, cond_init]; %integer(l), cond_lat, l>=1,
% cond = [integer(t), t>=0, t<=s_max, mem_all<=mem_tar, cond_init]; %, t(t_zero) <= 1
diagnostics = optimize(cond, obj, ops); %
disp(diagnostics);

if diagnostics.problem == 3
    % Maximum iterations exceeded
%     str = yalmiperror(diagnostics.problem);
%     warning(str);
elseif diagnostics.problem == 1
    % infeasible problem
    str = yalmiperror(diagnostics.problem);
    error(str);
elseif diagnostics.problem ~= 0
    % other errors
    str = yalmiperror(diagnostics.problem);
    error(str);
else
    % successflly solved
%     str = yalmiperror(diagnostics.problem);
%     disp(str);
end

% switch back when BNB does not produce better solution
if (cost_GP < value(cost) && (mem_GP <= mem_tar || ops.bnb.maxiter == 0))
   assign(t, t_GP) 
end


%% record solution
% assign(t,ceil(value(t)));
ts = value(t);
tm = [value(mem_all_real),value(mem_b_real')];
tc = value(cost);
td = diagnostics.problem;

% view solution
disp('=== Max tile size');
disp(s_max');
disp('=== Tile size selection :');
disp(ts');
% disp('=== Lattice sizes :');
% disp(value(l'));
disp('=== Target memory size:');
disp(mem_tar);
disp('=== Mem size:');
disp(tm)
disp('=== Minimum tile access cost:');
disp(cost_min)
disp('=== tile access cost:');
disp(tc)

end

