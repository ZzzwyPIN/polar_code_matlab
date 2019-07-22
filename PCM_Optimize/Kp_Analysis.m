%%
clear
clc
load('Pe_N256_snr3.2_R5.mat');
N=256;
Ng=12;
K = N/2+Ng;
Kp = 2:2:64;

[Ptmp,~] = sort(P);

Set_Ap = K-Kp;
Set_Astand = K-Kp/2;

for i = 1:length(Kp)
    delta_Pr(i) = sum(Ptmp(1:Set_Astand(i)))-sum(Ptmp(1:Set_Ap(i)));
end

semilogy(Kp,delta_Pr,'k-*');
xlabel('Kp');
grid on

%%
clear
load('Pe_N512_snr3_R5.mat');
N=512;
Ng=12;
K = N/2+Ng;
Kp = 2:2:128;

[Ptmp,~] = sort(P);

Set_Ap = K-Kp;
Set_Astand = K-Kp/2;

for i = 1:length(Kp)
    delta_Pr(i) = sum(Ptmp(1:Set_Astand(i)))-sum(Ptmp(1:Set_Ap(i)));
end

semilogy(Kp,delta_Pr,'k-*');
xlabel('Kp');
grid on

%%
clear
clc
load('Pe_N1024_snr3_R5.mat');
N=1024;
Ng=12;
K = N/2+Ng;
Kp = 2:4:256;

Ptmp = sort(P);

Set_Ap = K-Kp;
Set_Astand = K-Kp/2;

for i = 1:length(Kp)
    delta_Pr(i) = sum(Ptmp(1:Set_Astand(i)))-sum(Ptmp(1:Set_Ap(i)));
end

semilogy(Kp,delta_Pr,'k-*');
xlabel('Kp');
grid on

%%
clear
clc
load('Pe_N256_snr3.2_R5.mat');
N=256;
Ng=12;
K = N/2+Ng;
Kp = 2:2:64;

[Ptmp,~] = sort(P);

Set_Ap = K-Kp;
% Set_Astand = K-Kp/2;

for i = 1:length(Kp)-1
    delta_Pr(i) = sum(Ptmp(1:Set_Ap(i+1)))-sum(Ptmp(1:Set_Ap(i)));
end

semilogy(Kp,delta_Pr,'k-*');
xlabel('Kp');
grid on