
% nonlinear constraint for mpc

function [x,xeq] = xnonlcon(v) %% since only one input 'u', I need to write all informations inside again
global A;
global Beta;
global N;
global B_total;
global u_ot;
global x0;
global num_u;
global Num_region
global xeb
global num_ot

for t = 1:N
    x(Num_region*(t-1)+1:Num_region*t) =(B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)])-ones(Num_region,1); %% B*u < 1;
end

x1=x0;

for t = 1:N
    x1 = A*x1+Beta-(B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x1;
end

x(Num_region*N+1:Num_region*(N+1))= xeb-0.5.*ones(Num_region,1)-x1; %%terminal

x(Num_region*(N+1)+1:Num_region*(N+2))= -xeb-0.5.*ones(Num_region,1)+x1;
% x1
% xeb
% xeq =xeb-x1;
xeq = [];
% xeq(1:4)=v(N*num_u+1:N*num_u+Num_region,1)-x0;
% xeq(5:8)=A*v(N*num_u+1:N*num_u+Num_region,1)+Beta-(B*v(1:8,:)).*v(N*num_u+1:N*num_u+Num_region,1)-v(N*num_u+Num_region+1:N*num_u+2*Num_region,1);
%     for t1 = 1:N
%     xeq(t1*Num_region+1:t1*Num_region+Num_region)=A*v(N*num_u+(t1-1)*Num_region+1:N*num_u+t1*Num_region,1)+Beta-(B*v((t1-1)*num_u+1:t1*num_u,:)).*v(N*num_u+(t1-1)*Num_region+1:N*num_u+t1*Num_region,1)-v(N*num_u+t1*Num_region+1:N*num_u+t1*Num_region+Num_region,1);
%     end
end

