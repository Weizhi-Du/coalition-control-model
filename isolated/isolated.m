clc;
clear;
%% 初始化需要用到的变量——数据的预处理过程
% base parameters 基础参数
Num_region = 4; % number of region 区域数量
Num_boat = 6; % total number of boat; 渔船数量

%% region's parameters (no fishing)
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

% Fish inflow
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
%  x(:,1) = [1000, 1000, 1000, 1000]';
% no fishing dynamic
x_no(:,1) = x(:,1);
for t = 1:100
%     x_no(:,t+1) = A*x_no(:,t)+Beta; % xi(t+1)=Ai+Bixi(t)
        x_no(:,t+1) = A+Beta*x_no(:,t);
end

% boat's parameters
% gamma(1) = 0.08; 此处注释掉的为旧代码
% gamma(2) = 0.1;
% gamma(3) = 0.12;
% gamma(4) = 0.15;
% gamma(5) = 0.20;
% gamma(6) = 0.28;
% Gam = cell(6,1);
% Gam{1} = eye(4)*gamma(1);
% Gam{2} = eye(4)*gamma(2);
% Gam{3} = eye(4)*gamma(3);
% Gam{4} = eye(4)*gamma(4);
% Gam{5} = eye(4)*gamma(5);
% Gam{6} = eye(4)*gamma(6);
% 
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
N = 15; %prediction horizon 预测范围

%% isolated coalition 孤立联盟
global xeb;
max_iter = 720; %天数：720天，即迭代720次，对应paper中的time step
% 迭代求解，目的是得到渔船每天实际工作时的u值和其对应的捕鱼量
for t_d = 1:1:max_iter
    tic;
    x0 = x(:,end);
    %% 渔船1
    % 只针对渔船1，因为B=[Gam{1}],得到的B值，只和渔船1有关
    B = [Gam{1}];
    Num_alleffort = size(B,2);
    % 利用均衡函数去求渔区的均衡鱼量（xeb）和渔船的计划u值（veb）
    [veb,Jeq(1)] = equilibrium(A,Beta,B,Num_alleffort,Num_region,Num_boat,x0);
    xeb = veb(end-Num_region+1:end,1);
    % 将均衡函数求得的计划u值（effort，上面的veb）和xeb值（上面的均衡鱼量）带入mpc，得到优化的u值（渔船实际进行工作时的u值）；
    [v,J_sum,J]  = mpc_solver(A,Beta,N,B,x0,Num_region,xeb,veb);
    u1 = v(1:Num_alleffort,:); 
    
    %% 渔船2
    B = [Gam{2}];
    Num_alleffort = size(B,2);
    [veb,Jeq(2)] = equilibrium(A,Beta,B,Num_alleffort,Num_region,Num_boat,x0);
    xeb = veb(end-Num_region+1:end,1);
    [v,J_sum,J]  = mpc_solver(A,Beta,N,B,x0,Num_region,xeb,veb);
    u2 = v(1:Num_alleffort,:);
    %% 渔船3
    B = [Gam{3}];
    Num_alleffort = size(B,2);
    [veb,Jeq(3)] = equilibrium(A,Beta,B,Num_alleffort,Num_region,Num_boat,x0);
    xeb = veb(end-Num_region+1:end,1);
    [v,J_sum,J]  = mpc_solver(A,Beta,N,B,x0,Num_region,xeb,veb);
    u3 = v(1:Num_alleffort,:);
    %% 渔船4
    B = [Gam{4}];
    Num_alleffort = size(B,2);
    [veb,Jeq(4)] = equilibrium(A,Beta,B,Num_alleffort,Num_region,Num_boat,x0);
    xeb = veb(end-Num_region+1:end,1);
    [v,J_sum,J]  = mpc_solver(A,Beta,N,B,x0,Num_region,xeb,veb);
    u4 = v(1:Num_alleffort,:);
    %% 渔船5
    B = [Gam{5}];
    Num_alleffort = size(B,2);
    [veb,Jeq(5)] = equilibrium(A,Beta,B,Num_alleffort,Num_region,Num_boat,x0);
    xeb = veb(end-Num_region+1:end,1);
    [v,J_sum,J]  = mpc_solver(A,Beta,N,B,x0,Num_region,xeb,veb);
    u5 = v(1:Num_alleffort,:);   
    %% 渔船6
    B = [Gam{6}];
    Num_alleffort = size(B,2);
    [veb,Jeq(6)] = equilibrium(A,Beta,B,Num_alleffort,Num_region,Num_boat,x0);
    xeb = veb(end-Num_region+1:end,1);
    [v,J_sum,J]  = mpc_solver(A,Beta,N,B,x0,Num_region,xeb,veb);
    u6 = v(1:Num_alleffort,:);
    
    %将6艘船的u值拼成一个u值向量，以便进行后面的运算
    u = [u1;u2;u3;u4;u5;u6];
    B_total = [Gam{1} Gam{2} Gam{3} Gam{4} Gam{5} Gam{6}];
    for a_n = 1:1:Num_boat
    % 利用u值再求每艘船每天的捕鱼量，存入J_agent中
        J_agent(a_n,t_d) = ones(1,Num_region)*(B_total(1:Num_region,(a_n-1)*Num_region+1:a_n*Num_region)*u((a_n-1)*Num_region+1:a_n*Num_region,1).*x(:,end));
    end
    % 最后更新x，x是每天的实际鱼量，在下一轮求解时要用x0 = x(:,end);
    x(:,1+t_d) = Beta*x(:,end)+A-B_total*u.*x(:,end);
    toc;
end
%% function below
