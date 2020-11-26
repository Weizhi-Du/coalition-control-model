
function f = cost_hat(v,v_all,B,B_all)
% Input
% v: a num_coalition *num_region* N vector, a coalition all effort for every
% region in N times
% v_all: a (N * num_region*num_boat) vector, all boat effort for every
% region in N times
% B: a num_region * num_region matrix, a boat parameter, \gamma
% B_all: a num_region *(num_region * num_boat) matrix, all boat parameter,
% \gamma
% Output
% f: a const(negative number), return the cost for the current coalition(v)
f=0;
f_all = 0;
global A;
global Beta;
global N;
global x0;
global Num_region
num_u = size(B,2);
num_all = size(B_all,2);
x = x0;
Sum = ones(1,Num_region);% Num_region is the number of region
for t = 1:N
    f = -Sum*((B*v(t*num_u-num_u+1:t*num_u,:)).*x)+f;
    f_all = -Sum*((B_all*v_all(t*num_all-num_all+1:t*num_all,:)).*x)+f_all;
    x = A*x+Beta-(B_all*v_all(t*num_all-num_all+1:t*num_all,:)).*x;
end
end

