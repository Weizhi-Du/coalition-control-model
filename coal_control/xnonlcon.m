% nonlinear constraint for mpc
function [x,xeq] = xnonlcon(A,Beta,N,B_total,u_ot,x0,num_u,Num_region,xeb,num_ot,v) %% since only one input 'u', I need to write all informations inside again
for t = 1:N
x(Num_region*(t-1)+1:Num_region*t) =(B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)])-ones(Num_region,1); %% B*u < 1;
end
x1=x0;
for t = 1:N
    x1 = Beta*x1+A-(B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x1;
end
x(Num_region*N+1:Num_region*(N+1))= xeb-0.5.*ones(Num_region,1)-x1; %%terminal
x(Num_region*(N+1)+1:Num_region*(N+2))= -xeb-0.5.*ones(Num_region,1)+x1;
xeq = [];
end
