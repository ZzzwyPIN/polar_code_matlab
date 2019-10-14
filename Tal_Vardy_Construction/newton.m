function [x,iter]=newton(x0,f,fp) 
% newton-raphson algorithm 
MAX_RUN = 100;
N = MAX_RUN; eps = 1.e-5; % define max. no. iterations and error 
maxval = 10000.0; % define value for divergence 
xx = x0; 
while (N>0) 
    tmp1 = f(xx);
    tmp2 = fp(xx);
    if tmp2 == 0
        tmp1 = f(xx+0.01);
        tmp2 = fp(xx+0.01);
        xn = xx - tmp1 / tmp2;
    else
        xn = xx - tmp1  /tmp2;
    end
 if abs(f(xn))<eps 
     x=xn;iter=100-N; 
     return; 
 end; 
 if abs(f(xx))>maxval 
     disp(['iterations = ',num2str(iter)]); 
     error('Solution diverges'); 
     break; 
 end; 
 N = N - 1; 
 xx = xn; 
 iter=MAX_RUN-N;
end; 



