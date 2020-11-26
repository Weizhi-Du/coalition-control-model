%mpc solver

function [v,J_sum,J,u_to,u_toN,num_u] = mpc_slover(xeb,A,Beta,x0,Num_region,Gam,Num_boat,N,col,u_to,u_toN)
% Fixed other all coalition, Update the col current coalition
% N: a const, prediction horizon£¬max iteration num
% col: a num_coalition *1 vector, a coalition index vector
% u_to: a (num_boat * num_region)*1 vector, the all effort for boat
% u_toN: a (N * (num_boat * num_region)) * 1 vector, N all effort for boat
% Output
% v: a (N * num_coalition * num_region)*1 vector, all current coalition
% effort vector
% J_sum: a const, all current coalition fish for N times
% J: a const, fish carry num for one time current coalition fish
% u_to: Update (num_boat * num_region)*1 vector, update the all effort for boat
% u_toN:  Update (N * (num_boat * num_region)) * 1 vector,update N all effort for boat
% global x0;
% global Num_region;
% global B_total;
% global u_ot;
% global num_u;
% global Gam;
% global Num_boat;
% global num_ot
B = [];
c_else = [];
c_in = [];
v0 = [];
u_ot = [];%other coalition effort
for c = 1:1:size(col,2)
   if col(c) == 1  % loop the current coalition
       B = [B,Gam{c}];
       c_in = [c_in,c]; % the current coalition
   else 
       c_else = [c_else,c]; % the other coalition
   end
end

for vn = 1:1:N
    for cn = 1:1:size(c_in,2)% initial guess for all current coalition decision variable.
        v0 = [v0;u_toN(c_in(cn)*Num_region-Num_region+1+(vn-1)*Num_region*Num_boat:c_in(cn)*Num_region+(vn-1)*Num_region*Num_boat,:)]; 
    end
    for cn = 1:1:size(c_else,2)% initial guess for all current other coalition decision variable. 
        u_ot = [u_ot;u_toN(c_else(cn)*Num_region-Num_region+1+(vn-1)*Num_region*Num_boat:c_else(cn)*Num_region+(vn-1)*Num_region*Num_boat,:)];
    end
end

B_total = B; % give all boats B to B_total, then add other boates
if size(c_else,2) ~= 0    
    for ct = 1:1:size(c_else,2)
        B_total = [B_total,Gam{c_else(ct)}];
    end    
end

num_ru = size(B_total,2);%number of total efforts


num_u = size(B,2); %number of efforts in one time
num_ot = num_ru-num_u; % umber of other coalitions efforts.
num_b = num_u/Num_region; %number of boats
Sum = ones(1,Num_region);% Num_region is the number of region
% v0 = [kron(ones(N,1),u_to(1:num_u,:))];% initial guess for all decision variable. 
l_u = size(v0,1); %number of decision variable


% for c = 1:1:size(col,2)
%     if col(c) == 0
%     u_ot = [u_ot;u_to((c-1)*Num_region+1:c*Num_region,:)];
%     else
%     end
% end
% u_otN = kron(ones(N,1),u_ot);
    % lb = zeros(N*num_u,1);
    lb = zeros(l_u,1);
    ub =ones(N*num_u,1);

    Aeq = kron(ones(num_b*N,1),zeros(1,l_u)); % giving the size for equailty constraint.
    Aeq_u = kron(eye(num_b*N),ones(1,Num_region)); % giving the effort we need to calculate in equality constraint.
    Aeq(1:size(Aeq_u,1),1:size(Aeq_u,2))=Aeq_u;
    beq = ones(num_b*N,1);
%     nonlcon = @xnonlcon;
    options  = optimset( 'Algorithm','sqp','TolFun', 1e-5,'MaxIter',1000,'Tolcon',1e-5);
    v = fmincon(@(v)myfun(A,Beta,N,B_total,u_ot,x0,num_u,Num_region,num_ot,v),v0,[],[],Aeq,beq,lb,ub,@(v)xnonlcon(A,Beta,N,B_total,u_ot,x0,num_u,Num_region,xeb,num_ot,v),options);
    
    % give results to total effort(update); 
    v_st = 1; %initial an variable for counting how many boats in this current coalition
    for c = 1:1:size(col,2)
        if col(c) == 1
        u_to((c-1)*Num_region+1:c*Num_region,:) = v((v_st-1)*Num_region+1:v_st*Num_region,:);
        v_st = v_st+1;
        end
    end

    v_st = 1;
    for vn = 1:1:N
        for cn = 1:1:size(c_in,2)
            u_toN(c_in(cn)*Num_region-Num_region+1+(vn-1)*Num_region*Num_boat:c_in(cn)*Num_region+(vn-1)*Num_region*Num_boat,:)= v((v_st-1)*Num_region+1:v_st*Num_region,:); % update prediction
            v_st = v_st+1;
        end
    end
    
    J_sum = -myfun(A,Beta,N,B_total,u_ot,x0,num_u,Num_region,num_ot,v); % f: a const, all current coalition fish for N times
    J = Sum*((B*v(1:num_u)).*x0); % a const, one time current coalition fish
end

