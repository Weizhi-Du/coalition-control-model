function [v,Jeq] = equilibrium(A,Beta,B,Num_alleffort,Num_region,Num_boat,x0)
% function [v,Jeq] = equilibrium(x0)

% origin the effort
% B: a boat_num * boat_num matrix, B(i,i) is boat capability for i boat  
% v: origin the effort 
% Jeq: a const, O object function, the fish_num catch 
% global x0;
% global Num_region; 
% global Num_alleffort;  %number of efforts in one time
% Num_boat = Num_alleffort/Num_region; %number of boats  渔船数
Sum = ones(1,Num_region);% Num_region is the number of region
v0 = [0.25*ones(Num_alleffort,1);x0];% initial guess for all decision variable. 所有决策变量的初始猜测
l_u = size(v0,1); %number of decision variable  决策变量个数

fun = @(v)-Sum*((B*v(1:Num_alleffort,:)).*v(Num_alleffort+1:Num_region+Num_alleffort,:)); %objective function

% lb = zeros(N*Num_alleffort,1);
lb = zeros(l_u,1);
ub =ones(Num_alleffort,1);

Aeq = kron(ones(Num_boat,1),zeros(1,l_u)); % giving the size for equality constraint.
% Aeq_u = kron(eye(Num_boat),ones(1,Num_region)); % giving the effort we need to calculate in equality constraint.
% Aeq(1:size(Aeq_u,1),1:size(Aeq_u,2))=Aeq_u; 
beq = ones(Num_boat,1);
% nonlcon = @eqnonlcon;

options  = optimset( 'Algorithm','sqp','TolFun', 1e-5,'MaxIter',1000,'Tolcon',1e-5);
v = fmincon(fun,v0,[],[],Aeq,beq,lb,ub,@(v)eqnonlcon(A,Beta,B,Num_alleffort,Num_region,v),options); 
% v = fmincon(fun,v0,[],[],Aeq,beq,lb,ub,nonlcon,options);

Jeq = -Sum*((B*v(1:Num_alleffort,:)).*v(Num_alleffort+1:Num_region+Num_alleffort,1));
end

