function [x,xeq] = xnonlcon(A,Beta,N,B,x0,Num_alleffort,Num_region,xeb,v) %% since only one input 'u', I need to write all informations inside again
% restrict the fish number, to avoid destroying the fish balance
for t = 1:N
x(Num_region*(t-1)+1:Num_region*t) =(B*v(t*Num_alleffort-Num_alleffort+1:t*Num_alleffort,:))-ones(Num_region,1); %% B*u < 1;
end

x1=x0;
 
for t = 1:N
    x1 = Beta*x1+A-(B*v(t*Num_alleffort-Num_alleffort+1:t*Num_alleffort,:)).*x1;
end

x(Num_region*N+1:Num_region*(N+1))= xeb-0.5.*ones(Num_region,1)-x1; %%terminal
x(Num_region*(N+1)+1:Num_region*(N+2))= -xeb-0.5.*ones(Num_region,1)+x1;

xeq=[];

end

