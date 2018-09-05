%要证明规律3是否成立，只要证明最后一列的似然值是否已经解出。
clc;clear;
Np=9;
mem_lr = zeros(2^Np, Np+1);
for n=0:1:Np
    for j=1:1:2^(Np-n)
        mem_lr(2^n*(j-1)+1,Np+1-n) =1;
    end
end
mem_lr(2^(Np-1)+1,1)=1;

for j=2:2:(2^Np-2)
   idx2 = dec2bin(j);
   n_idx2 = length(idx2);
   for k=1:Np-n_idx2
       idx2 = strcat('0',idx2);
   end
   reverse_idx = bin2dec(fliplr(idx2))+1;
   node_idx = zeros(2^Np,Np+1);
   node_idx(1,1) = reverse_idx;
   node_idx(2,1) = reverse_idx+2^(Np-1);
   temp=0;
   for jj=2:1:Np+1 
        for jjj=1:1:2^(jj-2)
             %判断前后级节点之间的关系的方法，方法在ppt里已给出
             if mod(floor((node_idx(jjj,jj-1)-1)/2^(Np-jj+1)),2)==0
                  node_idx(2*jjj-1,jj) = node_idx(jjj,jj-1);
                  node_idx(2*jjj,jj) = node_idx(jjj,jj-1)+2^(Np-jj+1);
                  mem_lr(node_idx(2*jjj-1,jj),jj)=1;
                  mem_lr(node_idx(2*jjj,jj),jj)=1;
             else
                  temp=1;
                  node_idx(2*jjj-1,jj) = node_idx(jjj,jj-1)-2^(Np-jj+1);
                  node_idx(2*jjj,jj) = node_idx(jjj,jj-1);
             end
        end
        if temp==1 %若为左下蝶型则当前列已被解出，则跳出循环
            break
        end
   end
   k=jj;
   x=1;
   for jjj = 1:1:2^(k-1)
       if mem_lr(node_idx(jjj,k),k)==0
           x=0;
       end
   end
   if x==0
       break;
   end
end
if x==1
   disp('对于Np值，规律3是正确的。');
else
   disp('对于Np值，规律3是错误的。');
end