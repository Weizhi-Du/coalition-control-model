clc;
clear;
%%
global A;
global Beta;
global N;
global x0;
global Num_region
global Gam;
global Num_boat;
Num_region = 4; % number of region
Num_boat = 6; % total number of boat;

%% region's parameters (no fishing)
% Survival rate
% alpha(1) = 0.2;
% alpha(2) = 0.30;
% alpha(3) = 0.45;
% alpha(4) = 0.6;

r_alp = unifrnd(0.2, 0.6, Num_region);
for al = 1:Num_region
    alpha(al) = r_alp(1,al);
end

% Fish inflow
% beta(1) = 300;
% beta(2) = 450;
% beta(3) = 350;
% beta(4) = 200;
r_bet = randi([200 450],1,Num_region);
for be = 1:Num_region
    beta(be) = r_bet(1,be);
end

A = diag(alpha); % create a new diagonal matrix
Beta = beta';
% x(:,1) = [200, 300, 150, 250]'; % initial fish num
x(:,1) = randi([150 300],1,Num_region);

%% no fishing dynamic
x_no(:,1) = x(:,1);
for t = 1:100
    x_no(:,t+1) = A*x_no(:,t)+Beta; % xi(t+1)=Ai+Bixi(t)
end

%% boat's parameters
% gamma(1) = 0.08;
% gamma(2) = 0.1;
% gamma(3) = 0.12;
% gamma(4) = 0.15;
% gamma(5) = 0.20;
% gamma(6) = 0.28;

r_gamma = unifrnd(0.08, 0.28, Num_boat);
for ga = 1:Num_boat
    gamma(ga) = r_gamma(1,ga);
end
B_total=[];
Gam = cell(Num_boat,1);
for gam = 1:Num_boat
    Gam{gam} = eye(Num_region)*gamma(gam);
    B_total = [B_total,Gam{gam}];
end
% B_total = [Gam{1} Gam{2} Gam{3} Gam{4} Gam{5} Gam{6}];
% Gam{1} = eye(4)*gamma(1);
% Gam{2} = eye(4)*gamma(2);
% Gam{3} = eye(4)*gamma(3);
% Gam{4} = eye(4)*gamma(4);
% Gam{5} = eye(4)*gamma(5);
% Gam{6} = eye(4)*gamma(6);
% B = [Gam{1} Gam{2} Gam{3} Gam{4} Gam{5} Gam{6}];

% FORTEST max_it = 100;%max iteration for MPC;
max_it = 20;
% col_total = ones(1,Num_boat);% grand coalition
col_total =eye(Num_boat,Num_boat);

% col_total = diag(ones(Num_boat,1));
u_to = zeros(Num_region*Num_boat,1); %initial total effort

% FORTEST N = 15; %prediction horizon
N = 15;
Num_fish_Now = 1;
Num_fish_Last = 0;

%% cross communication for equilibrium
global xeb;
x0 = x(:,end);
xeb = x0;
max_iter_eq = 1; %max iteration algorithm for finding best equilibrium
% iterative algorithm
% For grand coalition, initialize the effort for every boat in different region
for t_d = 1:1:max_iter_eq
    xeb_c = [];
    for t_c = 1:1:size(col_total,1) % for every agent
        col = col_total(t_c,:);
        [veb,u_to] = equilibrium(col,u_to);
        xeb = veb(end-Num_region+1:end,1);
        xeb_c = [xeb_c;veb(end-Num_region+1:end,1)];
    end
    xeb_con(:,t_d) = xeb_c; % origin the steady state
    u_toc(:,t_d)=u_to;   % origin the effort vector for all boat
end
U_first =[];
U = [];
L = 1:Num_boat;
L = L';
for u_first = 1:Num_boat
    U_first = [U_first;(u_to(Num_region*(u_first-1)+1:(Num_region*u_first),1))'];  %将u拼成矩阵形式，每行代表一艘渔船的u值，每列为一个区域的u值
end
%% ACCM design
col_total = diag(ones(Num_boat,1)); %start with isolated coalition
col_new = {}; % for saving coalition formation
t_decide = 1; %decide to do operation of mereging ans splitting
u_toN = []; % for saving final MPC results in input signals
% FORTEST max_time = 720; %max time step of closed loop
max_time = 720;
u_cl = zeros(N*Num_region*Num_boat,max_time+10); % Save all medium effort status
fprintf('iterative begin，');
tic;
for t = 1:1:max_time
    x0 = x(:,end);
    disp(['iter = ',num2str(t)]);
    col_total = [];
    V_vir = [];
    % every 30 days to decide the merge/split
    if t_decide == 31 % decide merge/split after how many time step
        t_decide = 1;
    end
    if (t_decide == 1)  % get the current coalition status
        if size(unique(L),1) <= 3; %联盟数小于等于3时，不能再进行联盟
            col_total = col_new{t-1};
            [u_toN_in,vto_con_in,Jto_sum_in,J_sumc_in,J_hat_in] = convergempc(u_toc,col_total,x,xeb_con,max_it,V_vir,L);
            u_toN = u_toN_in;%output variable from spliting
            J_each_mpc(:,t) = J_hat_in(:,end);
            J_team_mpc(:,t) = J_sumc_in;
        else
            if t == 1
                U = U_first;
                V = U;
            else
                for pu = 1:Num_boat
                    U = [U;(u(Num_region*(pu-1)+1:(Num_region*pu),1))'];  %将u拼成矩阵形式，每行代表一艘渔船的u值，每列为一个区域的u值
                end
            end
            V_virtual = [];
            V_old = V; %此处的V值，是从上一轮迭代时求得的虚拟V值延续下来的，
                       %V值一直在变化，到这一轮时，先通过一个中间变量V_old把这个值记录下来，
                       %因为等下要和联盟以后的进行比较，如果联盟后的结果不好，还得采用这个旧的V值，即到时候再把V_old还给V.
            L_old = L; %L为每艘渔船的标签，如[1;1;3;3;5;5]表示1号船和2号船同属一个联盟，记为1（前一艘船的标号），3号船和4号船同属一个联盟3,5号船和6号船同属一个联盟5
                       %[1;1;1;1;5;5]表示1234号船合并成一个联盟，记为1，是由上一行中的{1;1}和{3;3}又进行了合并得到的，将{3;3}中的3号4号船标记成了1
            num = unique(L);  %unique for label 标签去重，得到的num向量为去重以后的标签向量
            k=size(num,1); %计算num的大小，得到的k即联盟个数
            col_total =zeros(k,Num_boat);
            for i=1:k
                w=find(L==num(i));%find the row num of ith col
                for j=1:size(w)
                    col_total(i,w(j))=1;
                end
            end
            [u_toN_in,vto_con_in,Jto_sum_in,J_sumc_in,J_hat_in] = convergempc(u_toc,col_total,x,xeb_con,max_it,V_vir,L);
            u_toN = u_toN_in;
            u = u_toN(1:Num_boat*Num_region);
            u_last = u;
            for a_n = 1:1:Num_boat
                J_agent_closed(a_n,t) = ones(1,Num_region)*(B_total(1:Num_region,(a_n-1)*Num_region+1:a_n*Num_region)*u((a_n-1)*Num_region+1:a_n*Num_region,1).*x(:,end));
            end
            J_agent_closed_Last = J_agent_closed;
            Num_fish_Last = sum(J_agent_closed(:,end));
           %% 合并操作准备：固定u求虚拟v
            m=1; %miu
            g=1; %gamma
            for o = 1:Num_boat % vk是一条船的，需要通过for循环把vk拼成最后的v
                V_o=find(L==L(o)); % 找到与第o个标签相同的船号(即行号)
                Pk = Num_boat - size(V_o,1); % Pk即与vl不在同一个联盟的船只数目
                V_SUM = zeros(1,Num_region);
                for h= 1:size(V_o,1)
                    % 计算与vo同一个联盟的船的v值之和
                    V_SUM = V_SUM + V(V_o(h),:);
                end
                V_SUM_PK = sum(V)-V_SUM; % 用全部船的v值之和减去刚计算的同一联盟的v值之和，即可得到与vo不在同一个联盟的船只的v值之和
                V_virtual_row = (1/(g*Pk-m))*(g*V_SUM_PK-m*U(o)); %计算第o艘船的虚拟v值
                V_virtual = [V_virtual;V_virtual_row]; % 拼成整体的虚拟v值
            end 
           %% 合并操作：利用虚拟v，进行合并
            [V,L] = itercol(V_virtual,L_old); %
            num = unique(L);
            k=size(num,1);
            col_total =zeros(k,Num_boat);
            for i=1:k
                w=find(L==num(i));
                for j=1:size(w)
                    col_total(i,w(j))=1;
                end
            end
            % 数据处理，没有实际物理含义，就是转换格式，变成后面的函数能用的格式
            for j = 1:size(V,1)
                TEMP = V(j,:);
                V_vir = [V_vir;TEMP'];
            end
            % 调用mpc
            [u_toN_in,vto_con_in,Jto_sum_in,J_sumc_in,J_hat_in] = convergempc(u_toc,col_total,x,xeb_con,max_it,V_vir,L);
            u_toN = u_toN_in;
            u = u_toN(1:Num_boat*Num_region);
            for a_n = 1:1:Num_boat
                J_agent_closed(a_n,t) = ones(1,Num_region)*(B_total(1:Num_region,(a_n-1)*Num_region+1:a_n*Num_region)*u((a_n-1)*Num_region+1:a_n*Num_region,1).*x(:,end));
            end
           %% 比较合并后的捕鱼量是否好于合并前的
            Num_fish_Now = sum(J_agent_closed(:,end));
            if Num_fish_Now < Num_fish_Last
                % 联盟后的效果不好，需要弃掉新的结果，继续沿用上一次迭代时的策略，即把刚才存的V_old、L_old、u_last和J_agent_closed_Last还给V、L、u和 J_agent_closed。
                V = V_old;
                L = L_old;
                u = u_last ;
                J_agent_closed =  J_agent_closed_Last;
            end % 如果 Num_fish_Now > Num_fish_Last，则无需进行额外处理，因为合并以后的u值和捕鱼量是后计算的，已经把旧的覆盖了
        end
    else % if(t_decide~=1)
        col_total = col_new{t-1};
        [u_toN_in,vto_con_in,Jto_sum_in,J_sumc_in,J_hat_in] = convergempc(u_toc,col_total,x,xeb_con,max_it,V_vir,L);
        u_toN = u_toN_in;%output variable from spliting
        J_each_mpc(:,t) = J_hat_in(:,end);
        J_team_mpc(:,t) = J_sumc_in;
    end
    t_decide = t_decide + 1;
    u = u_toN(1:Num_boat*Num_region);
    for a_n = 1:1:Num_boat
        J_agent_closed(a_n,t) = ones(1,Num_region)*(B_total(1:Num_region,(a_n-1)*Num_region+1:a_n*Num_region)*u((a_n-1)*Num_region+1:a_n*Num_region,1).*x(:,end));
    end
    Num_fish_Now = sum(J_agent_closed(:,end));
    x(:,1+t) = A*x(:,end)+Beta-B_total*u.*x(:,end);
    col_new{t} = col_total; % very import   
end
toc
Time = toc;
fprintf('iterative end, all time consumed(s)，');
disp(Time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%