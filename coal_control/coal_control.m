clc;
clear;
%% 初始化需要用到的变量——数据的预处理过程
% base parameters 基础参数
Num_region = 4; % number of region 区域数量
Num_boat = 6; % total number of boat; 渔船数量

%% region's parameters (no fishing)
% Survival rate
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
r_bet = randi([200 450],1,Num_region);
for be = 1:Num_region
    beta(be) = r_bet(1,be);
end

Beta = diag(alpha); % create a new diagonal matrix
A = beta';
% x(:,1) = [200, 300, 150, 250]'; % initial fish num
% x(:,1) = [1000, 1000, 1000, 1000]'; % initial fish num
x(:,1) = randi([150 300],1,Num_region);

%% no fishing dynamic
x_no(:,1) = x(:,1);
for t = 1:100
%     x_no(:,t+1) = A*x_no(:,t)+Beta; % xi(t+1)=Ai+Bixi(t)
        x_no(:,t+1) = A+Beta*x_no(:,t);
end

%% boat's parameters
% gamma(1) = 0.08;
% gamma(2) = 0.1;
% gamma(3) = 0.12;
% gamma(4) = 0.15;
% gamma(5) = 0.20;
% gamma(6) = 0.28;
% 初始化鱼流入量，定义每个区域的鱼流入量
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
max_it = 100;%max iteration for MPC;
% col_total = ones(1,Num_boat);% grand coalition
col_total =eye(Num_boat,Num_boat);

% col_total = diag(ones(Num_boat,1));
u_to = zeros(Num_region*Num_boat,1); %initial total effort
N = 15; %prediction horizon 

%% cross communication for equilibrium
% global xeb;
fprintf('equilibrium begin');
tic;
x0 = x(:,end);
xeb = x0;
max_iter_eq = 1; %max iteration algorithm for finding best equilibrium
% iterative algorithm
% For grand coalition, initialize the effort for every boat in different region
% 利用均衡函数去求渔区的均衡鱼量（xeb）和渔船的计划u值（veb）
for t_d = 1:1:max_iter_eq
    xeb_c = [];
    for t_c = 1:1:size(col_total,1) % for every agent/coalition
        col = col_total(t_c,:); 
        [veb,u_to] = equilibrium(A,Beta,xeb,Num_region,Gam,col,u_to);
        xeb = veb(end-Num_region+1:end,1);
        xeb_c = [xeb_c;veb(end-Num_region+1:end,1)];
    end
xeb_con(:,t_d) = xeb_c; % origin the steady state 
u_toc(:,t_d)=u_to;   % origin the effort vector for all boat
end
fprintf('equilibrium end，');
toc
Time_equilibrium = toc;
%% controlled coalition design 
col_total = diag(ones(Num_boat,1)); %start with isolated coalition
col_new = {}; % for saving coalition formation
t_decide = 1; %decide to do operation of mereging ans splitting
u_toN = []; % for saving final MPC results in input signals
max_time = 720; %max time step of closed loop
u_cl = zeros(N*Num_region*Num_boat,max_time+10); % Save all medium effort status
% num_u = size(B,2); %number of efforts in this coalition
tic;
fprintf('iterative begin，');
for t = 1:1:max_time 
    x0 = x(:,end);
    % every 30 days to decide the merge/split
    % 每30天考虑一次合并/分裂操作
    if t_decide == 31 % decide merge/split after how many time step
        t_decide = 1;    
    end
    if (t_decide == 1)  % get the current coalition status
    % merge operation
      if size(col_total,1)>1 % merge/split check when the decision day coming
            col_win = col_total;
            % compute the col_total effort, cost 
            [u_toN_in,vto_con_in,Jto_sum_in,J_sumc_in,J_hat_in] = convergempc(A,Beta,N,Num_region,Gam,Num_boat,u_toc,col_total,x,xeb_con,max_it);
            t_merge = 0;
            for col_1 = 1:1:size(col_total,1)-1 % for the col_1 coalition
                if col_1 > size(col_win,1) % if current index>all coalition index
                    break;
                end
                for col_2 = col_1+1:1:size(col_total,1) % for the col_2 coalition
                    if col_2 > size(col_win,1)%
                         break;
                    end
                    col_com = col_win;
                    col_com(col_1,:) = col_win(col_1,:)+col_win(col_2,:); % merge col_1 coalition and col_2 coalition
                    col_com(col_2,:)=[]; % delete the col_2 coalition
                    if (size(col_win,1) ~= size(col_total,1))&&(t_merge == 1)
                        % compute the isolated coalition effort
                        [u_toN_1,vto_con_1,Jto_sum_1,J_sumc_1,J_hat_1] = convergempc(A,Beta,N,Num_region,Gam,Num_boat,u_toc,col_total,x,xeb_con,max_it);
                        elseif size(col_win,1) == size(col_total,1)
                        u_toN_1 = u_toN_in;
                        vto_con_1 = vto_con_in;
                        Jto_sum_1 = Jto_sum_in;
                        J_sumc_1 = J_sumc_in;
                        J_hat_1 = J_hat_in;
                    end
                    % compute the merge effort 
                    [u_toN_2,vto_con_2,Jto_sum_2,J_sumc_2,J_hat_2] = convergempc(A,Beta,N,Num_region,Gam,Num_boat,u_toc,col_total,x,xeb_con,max_it);
                    t_merge = 0;
                    % the condition in selfish interest without redistribution            
                    Jm = 0;
                    Jn = 0;
                    JcomM = 0;
                    JcomN = 0;
                    for col_numk = 1:1:Num_boat
                        if col_win(col_1,col_numk)==1
                            Jm = J_hat_1(col_numk,end)+Jm;
                            JcomM = J_hat_2(col_numk,end)+ JcomM;
                        end
                        if col_win(col_2,col_numk)==1
                            Jn =  J_hat_1(col_numk,end)+Jn;
                            JcomN = J_hat_2(col_numk,end)+ JcomN;
                        end
                    end
                    % Statistics need to be integrated index 
                    if (JcomM>Jm) && (JcomN>Jn) % if the merge coalition> every isolated coalition
                        col_win = col_com;
                        t_merge = 1; % merge
                    end
                end
            end
           % final merge
           col_total = col_win;
           if t_merge == 1
                u_toN = u_toN_2;
           else
                u_toN = u_toN_1;
           end
        else
        [u_toN_in,vto_con_in,Jto_sum_in,J_sumc_in,J_hat_in] = convergempc(A,Beta,N,Num_region,Gam,Num_boat,u_toc,col_total,x,xeb_con,max_it);
        u_toN = u_toN_in;
      end
% splitting operation
    col_win = col_total;
    col_save_order = [];
    t_split = 0;
    for col_1 = 1:1:size(col_total,1)+Num_boat
        if col_1 > size(col_win,1)           
            break;         
        end     
        if (size(col_win,1) ~= size(col_total,1))&&(t_split == 1)
            % compute the split isoloted coalition
           [u_toN_1,vto_con_1,Jto_sum_1,J_sumc_1,J_hat_1] = convergempc(A,Beta,N,Num_region,Gam,Num_boat,u_toc,col_total,x,xeb_con,max_it);                              
        elseif size(col_win,1) == size(col_total,1) % if no split
           u_toN_1 = u_toN_in;
           vto_con_1 = vto_con_in;
           Jto_sum_1 = Jto_sum_in;
           J_sumc_1 = J_sumc_in;
           J_hat_1 = J_hat_in;    
        end
        d = sum(col_win(col_1,:));
        if d>1 % if including more than one boat
            [~,col] = find(col_win(col_1,:));
            for icol =2:1:size(col,2) %pick each boat
                col_save = col_win;
                col_save(col_1,col(1,icol))=0;
                col_save(size(col_win,1)+1,:)=zeros(1,size(col_win,2));
                col_save(size(col_win,1)+1,col(1,icol)) = 1; %save the boat               
                for i_save = 1:1:size(col_save,1)
                    pickcol= find(col_save(i_save,:),2);
                    col_save_order(i_save,1) = pickcol(:,1); 
                end
                colsort = sort(col_save_order); %ordering number pick       
                col_com = col_save; 
                % compute the split coalition index
                for i_save = 1:1:size(col_save,1)
                    rowold = find(col_save_order==colsort(i_save,:),1);                 
                    col_com(i_save,:) = col_save(rowold,:);% new coalition structure
                end
                col_2 = find(colsort==col_save_order(end,:),1);
                [u_toN_2,vto_con_2,Jto_sum_2,J_sumc_2,J_hat_2] = convergempc(A,Beta,N,Num_region,Gam,Num_boat,u_toc,col_total,x,xeb_con,max_it);
                t_split = 0; %reduce calculation when no change with the decision below:
                Jn = 0;
                JcomN = 0;
                for col_numk = 1:1:Num_boat
                    if col_com(col_2,col_numk)==1
                        Jn =  J_hat_1(col_numk,end)+Jn;
                        JcomN = J_hat_2(col_numk,end)+ JcomN;
                    end
                end
                if  JcomN> Jn
                    col_win = col_com;
                    t_split = 1;
                end
            end 
        end
    end
        if t_split == 1
            u_toN = u_toN_2;
            J_each_mpc(:,t) = J_hat_2(:,end);
            J_team_mpc(:,t) = J_sumc_2;
        else
            u_toN = u_toN_1;
            J_each_mpc(:,t) = J_hat_1(:,end);
            J_team_mpc(:,t) = J_sumc_1;
        end
        col_total = col_win;
    else % if(t_decide~=1)
    [u_toN_in,vto_con_in,Jto_sum_in,J_sumc_in,J_hat_in] = convergempc(A,Beta,N,Num_region,Gam,Num_boat,u_toc,col_total,x,xeb_con,max_it);
    u_toN = u_toN_in;%output variable from spliting
    J_each_mpc(:,t) = J_hat_in(:,end);
    J_team_mpc(:,t) = J_sumc_in;
    end
    t_decide = t_decide + 1;
    u = u_toN(1:Num_boat*Num_region);
    for a_n = 1:1:Num_boat
        J_agent_closed(a_n,t) = ones(1,Num_region)*(B_total(1:Num_region,(a_n-1)*Num_region+1:a_n*Num_region)*u((a_n-1)*Num_region+1:a_n*Num_region,1).*x(:,end));
    end
    x(:,1+t) = Beta*x(:,end)+A-B_total*u.*x(:,end);
    col_new{t} = col_total; % very import
end
fprintf('iterative end，');
toc
Time_iter = toc;
fprintf('all time consumed(s)，');
disp(Time_equilibrium+Time_iter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function below