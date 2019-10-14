function [y,lamda_next] = revcalcC( c,lamda0, sigma )

%options=optimoptions('fsolve','Algorithm',{'levenberg-marquardt',0.005},'Display','iter');
%x= fsolve(@(lamda)calcC_lamda(lamda,c),lamda0,options);

[x,iter] = newton(lamda0,@(lamda)calcC_lamda(lamda,c),@calcC_lamda_p);
%x is the lamda we required
%lamda = normpdf(y,1,sigma)/normpdf(y,-1,sigma) ==> lamda = exp(2y/sigma^2)
%y = x;
y = log(x)*sigma^2/2;
lamda_next = x;
end

