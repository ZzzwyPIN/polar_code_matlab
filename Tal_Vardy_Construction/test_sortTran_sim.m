% LL July 2, 2015
% test sortTran

L=4;

W=rand(1,2*L)

lr = W(1:L) ./ W(L+1:2*L)

Q = sortTran_sim(W)