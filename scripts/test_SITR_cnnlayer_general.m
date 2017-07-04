%% Symbolic Intra-Tile Memory Reuse

clear all;
bnb_maxiter = 300;

%% General CNN layer from FPGA'15 CNN paper (AlexNet)
% 

%% constants
% Layer 1: 317211
% input frame
N = 3;
inRow = 227;
inCol = 227;
% output frame
M = 48;
outRow = 55;
outCol = 55;
% convolution
K = 11; % kernel window size
ST = 4; % stride size

% % Layer 2 : 293040
% % input frame
% N = 48;
% inRow = 55;
% inCol = 55;
% % output frame
% M = 128;
% outRow = 27;
% outCol = 27;
% % convolution
% K = 5; % kernal window size
% ST = 1; % stride size

% % Layer 3: 532416
% % input frame
% N = 256;
% inRow = 16;
% inCol = 16;
% % output frame
% M = 192;
% outRow = 13;
% outCol = 13;
% % convolution
% K = 3; % kernal window size
% ST = 1; % stride size

% % Layer 4: 407424
% % input frame
% N = 192;
% inRow = 13;
% inCol = 13;
% % output frame
% M = 192;
% outRow = 13;
% outCol = 13;
% % convolution
% K = 3; % kernal window size
% ST = 1; % stride size

% % Layer 5: 286016
% % input frame
% N = 192;
% inRow = 13;
% inCol = 13;
% % output frame
% M = 128;
% outRow = 13;
% outCol = 13;
% % convolution
% K = 3; % kernal window size
% ST = 1; % stride size

%% loop iterator: (1)m, (2)r, (3)c, (4)n, (5)i, (6)j
n = 6;
s = sym('s', [n 1]);
s_min = 0*ones(n,1);
s_max = [M-1; outRow-1; outCol-1; N-1; K-1; K-1];
% s_min = [15;0;0;2;0;0];
% s_max = [15; outRow-1; outCol-1; 2; K-1; K-1];
assume(in(s, 'integer') & s>=s_min & s<=s_max)

%% iteration difference set/polytope: A*s <= b
A = zeros(2*n,n);
for i = 1:n
   A((2*i-1):2*i, i)  = [1;-1];
end
b = [reshape(s,1,n);reshape(s,1,n)];
b = reshape(b, 2*n, 1);
vtx = vtx_enum(A,b);

%% Array accesses
% in_layer access:
% [inCol*(n*inRow+(ST*r+i))+ST*c+j]
% [inCol*ST*r + ST*c + inCol*inRow*n + inCol*i + j]
% F_in = [0, inCol, 1, inCol*inRow, inCol, 1];
F_in = [0, 0, 0, 1, 0, 0;...
        0, ST, 0, 0, 1, 0;...
        0, 0, ST, 0, 0, 1];

% weight access
% [K*K*(m*N+n)+i*K+j)]
% [K*K*N*m + K*K*n + K*i + j]
F_wgt = [K*K*N, 0, 0, K*K, K, 1];

% out_layer access:
% [ourCol*(m*outRow+r)+c]
% [ourCol*outRow*r + ourCol*r + c]
F_out = [outCol*outRow, outCol, 1, 0, 0, 0];

%% reuse analysis
delete('log_cnnlayer_general_reuse.txt');
diary('log_cnnlayer_general_reuse.txt');
diary on;
disp('=================================================================');
disp('= In_layer access ');
disp('=================================================================');
rw_in = 1;
ker(1) = PRM( A, b, vtx, F_in, rw_in, s, s_min, s_max);

disp('=================================================================');
disp('= weight access ');
disp('=================================================================');
rw_wgt = 1;
ker(2) = PRM( A, b, vtx, F_wgt, rw_wgt, s, s_min, s_max);

disp('=================================================================');
disp('= Out_layer access ');
disp('=================================================================');
rw_out = 2;
ker(3) = PRM( A, b, vtx, F_out, rw_out, s, s_min, s_max);

diary off;

%% tile sizing with MIGP
% border: 16
diary on;
mem_tar = 2048; % 2048
[ ct, cost_rec, solu_rec, mem_rec, ker_rec ] = tile_size_opt( ker, A, b, vtx, n, s, s_min, s_max, mem_tar, bnb_maxiter );
diary off;

%% converge test
% diary on;
% nt = 11;
% mem_tar = [4096, 2^15];
% 
% for k=1:2
%     ct = zeros(nt,1);
%     cost_rec = zeros(nt,1);
%     solu_rec = zeros(n,nt);
%     mem_rec = zeros(nt,1);
% 
%     bnb_maxiter = 0;
%     for i = 1:nt
%         [ ct(i), cost_rec(i), solu_rec(:,i), mem_rec(i), ~ ] = tile_size_opt( ker, A, b, vtx, n, s, s_min, s_max, mem_tar(k), bnb_maxiter);
%         bnb_maxiter = bnb_maxiter + 30;
%     end
%     
%     str = sprintf('opt_cnnlayer3_converge_%d', mem_tar(k));
%     save(str, 'cost_rec', 'solu_rec', 'mem_rec', 'ct');
% end
% 
% diary off;


%% full tests
% diary on;
% nt = 13; % number of test cases
% 
% mem_tar = zeros(nt,1);
% ct = zeros(nt,1);
% cost_rec = zeros(nt,1);
% solu_rec = zeros(n,nt);
% mem_rec = zeros(nt,1);
% 
% for i = 1:nt
%     mem_tar(i) = 64*2^(i-1); % 2^6~18
%     [ ct(i), cost_rec(i), solu_rec(:,i), mem_rec(i), ~ ] = tile_size_opt( ker, A, b, vtx, n, s, s_min, s_max, mem_tar(i), bnb_maxiter );
% end
% 
% save('opt_cnnlayer1_mem', 'cost_rec', 'solu_rec', 'mem_rec', 'ct', 'mem_tar');
% 
% diary off;

%%
% diary on
% 
% disp('=================================================================');
% disp('= Full Optimiztion Test Result:');
% disp('=================================================================');
% 
% disp('=== solutions');
% disp('* max:');
% disp(s_max');
% disp('* records:');
% disp(solu_rec);
% 
% disp('=== mem size');
% disp('* target:');
% disp(mem_tar);
% disp('* records:');
% disp(mem_rec);
% 
% disp('=== tile access cost');
% disp(cost_rec);
% 
% disp('=== computation time:  ');
% disp(ct);
% 
% diary off



