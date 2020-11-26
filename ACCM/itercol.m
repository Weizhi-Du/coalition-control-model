function [ V,L ] = itercol( V,L )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%    L=[1;1;3;3;5;5];
num = unique(L);  %取去重以后的联盟标签
k=size(num,1);
Data=cell(2,k);

for i=1:k
     w=find(L==num(i));%找到第i个联盟所有船对应的船号（行号）
     A=[];
     Y=[];
     for j=1:size(w)
         [I]=V(w(j),:);%在原始船只矩阵中取出对应行号的船只行数据，放入A中
         [X]=w(j);
          A=[A;I];
          Y=[Y;X];
     end
     Data{1,i} = A;
     Data{2,i} = Y;
end

V_Temp = [];
for g = 1:k
    V_Temp = [V_Temp;Data{1,g}(1,:)];
end
D = pdist2(V_Temp,V_Temp);  %计算联盟间的距离
Dnew = D; 
Dnew(Dnew==0)=inf;
[minv,ind]=min(Dnew,[],2);
min_row = min(minv);
dmin = min_row;
% d=0.1;  %距离间隔设为0.1
while min_row == dmin
    z=find(minv==min_row); %发现距离最近的联盟号
    if size(z,1)<2
        break;
    else
        col_left = Data{1,z(1)}(1,:);
        col_right = Data{1,z(2)}(1,:);
        boatnum_left = size(Data{2,z(1)},1);
        boatnum_right = size(Data{2,z(2)},1);
        %计算合并后的虚拟v
        v_temp = (col_left*boatnum_left+col_right*boatnum_right)/(boatnum_left+boatnum_right);
        for left_i = 1:boatnum_left
            l1 = Data{2,z(1)}(left_i,:);
            V(l1,:) = v_temp; %替换成合并后的v
            L(l1,:) = Data{2,z(1)}(1,:); %更改联盟号
        end
        for right_i = 1:boatnum_right
            l2 = Data{2,z(2)}(right_i,:);
            V(l2,:) = v_temp;
            L(l2,:) = Data{2,z(1)}(1,:);
        end
        minv(z(1)) = inf;
        minv(z(2)) = inf;
        min_row = min(minv);
    end
end

end

