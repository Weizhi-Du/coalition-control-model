%iteration algorithm for mpc
function [u_toN,vto_con,Jto_sumc,J_sum,J_hat] = convergempc(A,Beta,N,Num_region,Gam,Num_boat,u_toc,col_total,x,xeb_con,max_it)
%
% Input
% u_toc: a (num_boat * num_region)*1 vector, initialize the all effort for boat
% col_total: a (num_boat * num_boat) matrix, coalition index matrix, could get the size
% x: a num_region * 1 vector, origin the fish num for every origin
% xeb_con: a num_region * 1 vector, compute the origin fish num
% max_it: a const, max iteration num
% Output
% u_toN: a (N* (num_boat * num_region))* 1 vector, N all effort vector
% vto_con: :a (num_coalition * N* (num_boat * num_region)) * num_boat, all coalition effort vector for N times
% Jto_sumc:a const, all coalition fish carry num for N times
% J_sum: a const, medium variable, no use 
% J_hat: a num_boat * num_boat matrix, No reallocation of resources, the
% number of fish acquired by each bot
% global N;
% global Num_region
% global Gam;
% global Num_boat;
% global x0;
% global num_u
u_to = u_toc(:,end);
u_toN = [kron(ones(N,1),u_to)];
num_ru = size(col_total,2)*Num_region; % coalition * Num_region
for t_d = 1:1:max_it
    v_con = [];
    J_sumc = [];
    x0 = x(:,end);
    for t_c = 1:1:size(col_total,1)
        col = col_total(t_c,:);
                 xeb = xeb_con((t_c-1)*Num_region+1:t_c*Num_region,end);
%         xeb = xeb_con(end-3:end,end);
        [v,J_sum,~,u_to,u_toN,num_u]  = mpc_slover(xeb,A,Beta,x0,Num_region,Gam,Num_boat,N,col,u_to,u_toN);
        u_out(:,t_c) = [v(1:num_u,:);zeros(num_ru-num_u,1)];
        v_con = [v_con;u_toN];
        J_sumc = [J_sumc;J_sum];
    end
    B_all = [];
    for B_c = 1:1:Num_boat
        B_all = [B_all,Gam{B_c}];
    end
    for t_hat = 1:1:Num_boat
        v_hat = [];
        for vn = 1:1:N
            v_hat =[v_hat;u_toN(t_hat*Num_region-Num_region+1+(vn-1)*Num_region*Num_boat:t_hat*Num_region+(vn-1)*Num_region*Num_boat,:)];
        end
        J_hat(t_hat,t_d) = -cost_hat(A,Beta,N,x0,Num_region,v_hat,u_toN,Gam{t_hat},B_all);
    end
    vto_con(:,t_d) = v_con;
    Jto_sumc(:,t_d) = J_sumc;
    if t_d > 2
        if sum(abs(J_hat(:,t_d)-J_hat(:,t_d-1))) < 1*10^-15
            break
        end
    end
    %     u = [u1;u2;u3;u4;u5;u6];
    %     B_total = [Gam{1} Gam{2} Gam{3} Gam{4} Gam{5} Gam{6}];
    %     x(:,1+t_d) = A*x(:,end)+Beta-B_total*u.*x(:,end);
end
end
