%equlibrium algorithm

function [v,u_to] = equilibrium(col,u_to)
% origin the effort for every boat in different region
% Input
% col: for col_th coalition
% u_to: a num_boat * num_region vector, initial total effort
% Output
% v: (a num_boat * num_region + num_region) vector ,the first num_boat *
%    num_region is the effort vector, the last num_region is fish num
% u_to: a num_boat * num_region vector, the update effort vector
global num_u;
global xeb;
global Num_region;
global B_total;
global u_ot;%other coalition effort
global Gam;

B = [];
c_else = [];
for c = 1:1:size(col,2)
    if col(c) == 1
        B = [B,Gam{c}];
    else
        c_else = [c_else,c];
    end
end

B_total = B; % give all boats B to B_total, then add other boates

if size(c_else,2) ~= 0
    for ct = 1:1:size(c_else,2)
        B_total = [B_total,Gam{c_else(ct)}];
    end
end

num_ru = size(B_total,2);%number of total efforts

num_u = size(B,2); %number of efforts in this coalition
num_b = num_u/Num_region; %number of boats
Sum = ones(1,Num_region);% Num_region is the number of region


u_ot = [];%other coalition effort
u_in = []; %this coalition previous effort
for c = 1:1:size(col,2)
    if col(c) == 0
        u_ot = [u_ot;u_to((c-1)*Num_region+1:c*Num_region,:)];
    else
        u_in=[u_in;u_to((c-1)*Num_region+1:c*Num_region,:)];
    end
end

if sum(u_in)==0
    v0 = [0.25*ones(num_u,1);xeb];% initial guess for all decision variable.
else
    v0 = [u_in;xeb];
end
l_u = size(v0,1); %number of decision variable


fun = @(v)-Sum*((B_total*[v(1:num_u,:);u_ot]).*v(num_u+1:Num_region+num_u,:)); %objective function

% lb = zeros(N*num_u,1);
lb = zeros(l_u,1);
ub =ones(num_u,1);

Aeq = kron(ones(num_b,1),zeros(1,l_u)); % giving the size for equailty constraint.
Aeq_u = kron(eye(num_b),ones(1,Num_region)); % giving the effort we need to calculate in equality constraint.
Aeq(1:size(Aeq_u,1),1:size(Aeq_u,2))=Aeq_u;
beq = ones(num_b,1);
nonlcon = @eqnonlcon;
options  = optimset( 'Algorithm','sqp','TolFun', 1e-5,'MaxIter',1000,'Tolcon',1e-5);

v = fmincon(fun,v0,[],[],Aeq,beq,lb,ub,nonlcon,options);


% give results to total effort(update);

v_st = 1; %initial an variable for counting how many boats in this coalition
for c = 1:1:size(col,2)
    if col(c) == 1
        u_to((c-1)*Num_region+1:c*Num_region,:) = v((v_st-1)*Num_region+1:v_st*Num_region,:);
        v_st = v_st+1;
    end
end
end

