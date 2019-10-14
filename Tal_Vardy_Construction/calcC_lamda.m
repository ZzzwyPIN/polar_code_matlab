function c_lamda = calcC_lamda( lamda,c)
% C as in (33)
%   
%lamda as function of y
%I use 1 and -1 as centers to comply with author's requirement
%lamda = normpdf(y,1,sigma)/normpdf(y,-1,sigma);
% now calculate C according to formula in the paper
c_lamda = 1 - lamda* log2(1+1/lamda)/(1+lamda) - log2(lamda+1)/(1+lamda)-c;
end

