function [v,J_sum,J] = mpc_solver(A,Beta,N,B,x0,Num_region,xeb,veb)
% Solve the optimal schedule
%Input:
% N: a const, prediction horizon£¬max iteration num
% B: a boat_num * boat_num matrix, B(i,i) is i_boat capability
% veb: a (Num_alleffort+ Num_fish_region) * 1 vector, medium variable
%Output:
% v: 28*1,v(1:Num_alleffort24) is the current effort
% J_sum: a const, J_sum is the num of all region (all time,15 steps)
% J: is the fish catches (one time)
Num_alleffort = size(B,2); %number of efforts in one time
Num_boat = Num_alleffort/Num_region; %number of boats
v0 = [kron(ones(N,1),veb(1:Num_alleffort,1))];% initial guess for all decision variable.
l_u = size(v0,1); %number of decision variable
lb = zeros(l_u,1); %lower bound
ub = ones(l_u,1);
% compute the matrix to guarentee the sum_effort is 1
Aeq = kron(ones(Num_boat*N,1),zeros(1,l_u)); % giving the size for equality constraint.
Aeq_u = kron(eye(Num_boat*N),ones(1,Num_region)); % giving the effort we need to calculate in equality constraint.
Aeq(1:size(Aeq_u,1),1:size(Aeq_u,2))=Aeq_u;
beq = ones(Num_boat*N,1);
options = optimset( 'Algorithm','sqp','TolFun', 1e-5,'MaxIter',1000,'Tolcon',1e-5);
v = fmincon(@(v)myfun(A,Beta,N,B,x0,Num_alleffort,Num_region,v),v0,[],[],Aeq,beq,lb,ub,@(v)xnonlcon(A,Beta,N,B,x0,Num_alleffort,Num_region,xeb,v),options);
J_sum = -myfun(A,Beta,N,B,x0,Num_alleffort,Num_region,v);
J = sum(((B*v(1:Num_alleffort)).*x0));
end