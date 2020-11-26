%nonlinear constraint for equilibrium algorithm

function [x,xeq] = eqnonlcon(v) %% since only one input 'u', I need to write all informations inside again
global A;
global Beta;
global B_total;
global num_u;
global Num_region
global u_ot;

x =[];
xeq(1:Num_region)=A*v(num_u+1:Num_region+num_u,:)+Beta-(B_total*[v(1:num_u,:);u_ot]).*v(num_u+1:Num_region+num_u,:)-v(num_u+1:Num_region+num_u,:);

end

