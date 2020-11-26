%objective function for mpc algorithm

function f = myfun(A,Beta,N,B_total,u_ot,x0,num_u,Num_region,num_ot,v)
% f: a const, all current coalition fish for N times
f=0;
% 
% global A;
% global Beta;
% global N;
% global B_total;
% global u_ot;
% global x0;
% global num_u;
% global Num_region
% global num_ot

x = x0;
Sum = ones(1,Num_region);% Num_region is the number of region
    for t = 1:N
        f = -Sum*((B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x)+f;
        x = Beta*x+A-(B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x;
    end
end

