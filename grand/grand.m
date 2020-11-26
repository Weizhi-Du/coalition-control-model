clc;
clear;
%% 初始化需要用到的变量——数据的预处理过程
% base parameters 基础参数
Num_region = 4; % number of region 区域数量
Num_boat = 6; % total number of boat; 渔船数量

% region's parameters
% Survival rate 此处注释掉的为旧代码
% alpha(1) = 0.2;
% alpha(2) = 0.30;
% alpha(3) = 0.45;
% alpha(4) = 0.6;

% 初始化变量存活率，定义鱼在每个区域的存活率
r_alp = unifrnd(0.2, 0.6, Num_region);
%随机生成0.2-0.6之间的数，生成一个全是随机数的矩阵，命名为r_alp
% unifrnd 和 底下的randi是matlab封装好的函数，直接当黑盒子用，知道其是为了生成随机数即可
for al = 1:Num_region
    alpha(al) = r_alp(1,al); %遍历（循环访问）区域，给每个区域依次赋值，值为r_alp矩阵中的元素
end

% Fish inflow 此处注释掉的为旧代码
% beta(1) = 300;
% beta(2) = 450;
% beta(3) = 350;
% beta(4) = 200;
% 初始化鱼流入量，定义每个区域的鱼流入量
r_bet = randi([200 450],1,Num_region); 
for be = 1:Num_region
    beta(be) = r_bet(1,be);
end

Beta = diag(alpha); % create a new diagonal matrix
A = beta';
% x(:,1) = [200, 300, 150, 250]'; % initial fish num
x(:,1) = randi([150 300],1,Num_region);
% x(:,1) = [1000, 1000, 1000, 1000]';
% no fishing dynamic
x_no(:,1) = x(:,1);
for t = 1:100
%     x_no(:,t+1) = A*x_no(:,t)+Beta; % xi(t+1)=Ai+Bixi(t)
        x_no(:,t+1) = A+Beta*x_no(:,t);
end

% boat's parameters 此处注释掉的为旧代码
% gamma(1) = 0.08;
% gamma(2) = 0.1;
% gamma(3) = 0.12;
% gamma(4) = 0.15;
% gamma(5) = 0.20;
% gamma(6) = 0.28;
% r_gamma = unifrnd(0.08, 0.28, Num_boat);
r_gamma = unifrnd(0.01, 0.05, Num_boat);
for ga = 1:Num_boat
    gamma(ga) = r_gamma(1,ga);
end
B_total=[];
Gam = cell(Num_boat,1);
for gam = 1:Num_boat
    Gam{gam} = eye(Num_region)*gamma(gam);
    B_total = [B_total,Gam{gam}];
end
N = 15; %prediction horizon

%% grand coalition 大联盟
%% 均衡函数
% global xeb;
fprintf('grand coalition begin');
fprintf('\n');
B = B_total;
x0 = x(:,1);
Num_alleffort = Num_boat * Num_region; %number of efforts in one time
fprintf('equilibrium begin');
tic;
% <--equilibrium，using fmincon（a fun，a eqnonlcon，% simple compute the initial value ）-->
% 利用均衡函数去求渔区的均衡鱼量（xeb）和渔船的计划u值（veb）
[veb,Jeq]= equilibrium(A,Beta,B,Num_alleffort,Num_region,Num_boat,x0); 
xeb = veb(Num_alleffort+1:end,1);
% <------------------------------------------------------------------->
fprintf('equilibrium end，');
toc
Time_equilibrium = toc;

%% 迭代求解mpc
fprintf('iterative mpc_solver');
tic;
max_iter = 720; %天数：720天，即迭代720次，对应paper中的time step
% 迭代求解，目的是得到渔船每天实际工作时的u值和其对应的捕鱼量
for t_d = 1:1:max_iter % time
    tic
    x0 = x(:,end);
    % 将均衡函数求得的计划u值（effort，上面的veb）和xeb值（上面的均衡鱼量）带入mpc，得到优化的u值（渔船实际进行工作时的u值）；
    [v,J_sum(t_d),J(t_d)]  = mpc_solver(A,Beta,N,B,x0,Num_region,xeb,veb); % J is the fish catches
    u = v(1:Num_alleffort,1); % effort
    for a_n = 1:1:Num_boat % number of boats
    % 利用u值再求每艘船每天的捕鱼量，存入J_agent中
        J_agent(a_n,t_d) = sum(B_total(1:Num_region,(a_n-1)*Num_region+1:a_n*Num_region)*u((a_n-1)*Num_region+1:a_n*Num_region,1).*x(:,end));
    end
    % 最后更新x，x是每天的实际鱼量，在下一轮求解时要用x0 = x(:,end);
    x(:,1+t_d) = Beta*x(:,end)+A-B*u.*x(:,end);
    toc
end
fprintf('iterative end，');
toc
Time_mpc_solver = toc;
fprintf('grand coalition end，time(s)：');
Time_grand = Time_equilibrium + Time_mpc_solver;
disp(Time_grand);
