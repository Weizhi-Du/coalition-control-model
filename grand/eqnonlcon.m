function [x,xeq] = eqnonlcon(A,Beta,B,Num_alleffort,Num_region,v) %% since only one input 'u', I need to write all informations inside again
x = [];
xeq(1:Num_region)=Beta*v(Num_alleffort+1:Num_region+Num_alleffort,:)+A-(B*v(1:Num_alleffort,:)).*v(Num_alleffort+1:Num_region+Num_alleffort,:)-v(Num_alleffort+1:Num_region+Num_alleffort,:);
end

