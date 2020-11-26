 function [x,xeq] = eqnonlcon(A,Beta,B,Num_alleffort,Num_region,v) %% since only one input 'u', I need to write all informations inside again
% function [x,xeq] = eqnonlcon(v) %% since only one input 'u', I need to write all informations inside again
% global A;
% global Beta;
% global B;
% global Num_alleffort;
% global Num_region

x = [];
% xeq(1:Num_region)=A*v(Num_alleffort+1:Num_region+Num_alleffort,:)+Beta-(B*v(1:Num_alleffort,:)).*v(Num_alleffort+1:Num_region+Num_alleffort,:)-v(Num_alleffort+1:Num_region+Num_alleffort,:);
xeq(1:Num_region)=Beta*v(Num_alleffort+1:Num_region+Num_alleffort,:)+A-(B*v(1:Num_alleffort,:)).*v(Num_alleffort+1:Num_region+Num_alleffort,:)-v(Num_alleffort+1:Num_region+Num_alleffort,:);

end

