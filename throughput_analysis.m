% LL May 27, 2019
% justify the throughput of PCM
% 应该用throughput/area来作为衡量
clc; clear; close all;

Rs = 100; % in Mbps, source generation
R = 1/2; % the code rate
N=256; % the code length
K = floor(N*R); % information bits in one block
Kp = 32; % mutual information bits among two adjecent blocks
timeScale = 8;
Ts = 1/Rs; % the symbol duration in us

alpha = 0.1:1e-1:1; % decoding time
Td = alpha * Ts; % decoding time in proportional of the symbol duration
% Td = alpha * Ts * N;
% variable to control the second round decoding time
% can be adjusted
beta = 0.42;
Td_additional = Td * beta; % 第一张图上数值多少就是多少
% frozenbit_numPCM = 2*N-2*K+Kp;


% 2.5 dB
PB = 2e-3;
% throughput analysis
% case 1
% R1 = (1-PB)^2*(2*K-Kp+frozenbit_numPCM) ./ (2*N*Ts+2*Td); %  冻结位应该也算进去。刘强
R1 = (1-PB)^2*2*N./(2*N*Ts+2*Td);
% R1 = (1-PB)^2*2*N./(2*Td);
% case 2/3
% R23 = 2*(1-PB)*PB*(2*K-Kp+frozenbit_numPCM)./ (2*N*Ts+2*Td+Td_additional);
R23 = 2*(1-PB)*PB*2*N./(2*N*Ts+2*Td+Td_additional);
% R23 = 2*(1-PB)*PB*2*N./(2*Td+Td_additional);

R_new = R1+R23; 

R_new_norm = R_new / Rs;


% SCL, L=8
PB_scl = 1e-3;
% the decoding time in terms of the SC decoding
Td_scl = Td * timeScale; %1.5
R_scl = 0.4375;
K_scl = floor(N*R_scl);
% frozenbit_numSCL = N-K_scl;
% R_new_scl = ((1-PB_scl)*K_scl+frozenbit_numSCL)./(N*Ts + Td_scl);
R_new_scl = (1-PB_scl)*N./(N*Ts + Td_scl);
% R_new_scl = (1-PB_scl)*N./(Td_scl);
R_scl_norm = R_new_scl / Rs;
plot(alpha, R_new_norm,'k-o', alpha, R_scl_norm, 'r-*');
legend('PCM', 'SCL8(Tscl = 8Tsc)');
xlabel('T_{sc}/T_s');
ylabel('Normalized Throughput');
% axis([.1 1 .98 1]);
name = ['T_{scl} = ',num2str(timeScale),'T_{sc}'];
title(name)