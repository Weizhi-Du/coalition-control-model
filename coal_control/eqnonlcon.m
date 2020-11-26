%nonlinear constraint for equilibrium algorithm

function [x,xeq] = eqnonlcon(A,Beta,B_total,num_u,Num_region,u_ot,v) %% since only one input 'u', I need to write all informations inside again
% global A;
% global Beta;
% global B_total;
% global num_u;
% global Num_region
% global u_ot;

x =[];
xeq(1:Num_region)=Beta*v(num_u+1:Num_region+num_u,:)+A-(B_total*[v(1:num_u,:);u_ot]).*v(num_u+1:Num_region+num_u,:)-v(num_u+1:Num_region+num_u,:);

end

