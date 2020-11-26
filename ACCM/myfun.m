%objective function for mpc algorithm

function f = myfun(V_vir,L,v)
% f: a const, all current coalition fish for N times
f=0;

global A;
global Beta;
global N;
global B_total;
global u_ot;
global x0;
global num_u;
global Num_region
global num_ot

x = x0;
Sum = ones(1,Num_region);% Num_region is the number of region
% 旧myfun函数的公式。变量V_vir为NULL时启用，即没有进行合并操作
if isempty(V_vir)
    for t = 1:N
    f = -Sum*((B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x)+f;
    x = A*x+Beta-(B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x;
    end
else
   % 新myfun函数的公式。求出虚拟V_vir时启用，即进行了合并操作
   m = 0.001; % miu值
   z = unique(L);    
   for t = 1:length(z)
       LIE = find(L==z(t));
       vir_TEMP = [];
       for f = 1:size(LIE)
           h = LIE(f);
           vir_TEMP = [vir_TEMP;V_vir(h*Num_region-Num_region+1:h*Num_region,1)];
       end
       V_vir_NEW = repmat(vir_TEMP,15,1);
       f = -Sum*((B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x)+m*sum((v(t*num_u-num_u+1:t*num_u,:)-V_vir_NEW(t*num_u-num_u+1:t*num_u,:)).^2)+f;
       %         f = -Sum*((B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x)+f;
       x = A*x+Beta-(B_total*[v(t*num_u-num_u+1:t*num_u,:);u_ot(t*num_ot-num_ot+1:t*num_ot,:)]).*x;
   end
end
end

