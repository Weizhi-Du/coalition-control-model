function f = myfun(A,Beta,N,B,x0,Num_alleffort,Num_region,v)
% object function, loss function, current total fish caught number for all region
f=0;
x = x0;
Sum = ones(1,Num_region);% Num_region is the number of region
for t = 1:N
    f = -Sum*((B*v(t*Num_alleffort-Num_alleffort+1:t*Num_alleffort,:)).*x)+f; % total fish caught
    x = Beta*x+A-(B*v(t*Num_alleffort-Num_alleffort+1:t*Num_alleffort,:)).*x; % number of fish
end
end


