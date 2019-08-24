

%%
clear
clc
load('PCM_K140_CRC12_Kp24_R4531.mat');


figure
semilogy(SNR,berSC,'b-+',SNR,berBP,'k-*',SNR_SCL,berSCL2,'r-s',SNR_SCL,berSCL4,'m-^');
hold on
semilogy(SNR_PCM(1:4),berPCM_SC,'k-o',SNR_PCM(1:4),berPCM_BP,'b-p');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('BER','interpreter','Latex');
axis([1 4.5 1.0e-5 1])
legend('SC','BP','SCL, $L=2$','SCL, $L=4$','PCM-SC-2','PCM-BP-2','interpreter','Latex');
grid on

% alpha_b = 0.38;
% alpha_c = 6.9;
% perNewMax = (1+alpha_b).*P_B.^2 - alpha_b.*P_B.^3;
% perNewMin = (1+alpha_c).*P_B.^2 - alpha_c.*P_B.^3;

figure
semilogy(SNR,perSC,'b-+',SNR,perBP,'k-*',SNR_SCL,perSCL2,'r-s',SNR_SCL,perSCL4,'m-^');%
hold on
semilogy(SNR_PCM,perPCM_SC,'k-o',SNR_PCM,perPCM_BP,'b-p');%
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('PER','interpreter','Latex');
axis([1 4.5 1.0e-4 1])
legend('SC','BP','SCL, $L=2$','SCL, $L=4$','PCM-SC-2','PCM-BP-2','interpreter','Latex');
grid on

%% 
clear
clc
load('PCM_K140_CRC12_Kp32_R4375.mat');
figure
semilogy(SNR,perSC,'k-+',SNR(1:5),perBound,'k-h');
hold on
semilogy(SNR(1:5),perPCM_SC_2,'r-o');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('Packet Error Rate');
% axis([2 5 1.0e-5 1])
legend('SC','perBound','PCM-SC-2');
grid on

%%
clear
clc
load('PCM_SCL_K140_Kp32_CRC12_R4375.mat')
figure
semilogy(SNR,perSCL4,'b-+',SNR,perSCL8,'c-+',SNR,perSCL16,'m-h');
hold on
semilogy(SNR,per2,'k-^',SNR,per4,'r-o',SNR(1:4),per8,'b-s');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('Packet Error Rate');
% axis([2 5 1.0e-5 1])
legend('SCL,L=4','SCL,L=8','SCL,L=16','PCM-SCL,L=2','PCM-SCL,L=4','PCM-SCL,L=8');
grid on


%%
clear
clc
load('PCM_SCL_K140_Kp24_CRC12_R4531.mat')
figure
semilogy(SNR,perSCL4,'b-+',SNR,perSCL8,'c-^',SNR2,perSCL16,'m-h');
hold on
semilogy(SNR,per2,'k-*',SNR2,per4,'r-o',SNR2,per8,'b-s');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('PER','interpreter','Latex');
axis([1 3.5 1.0e-5 1])
legend('SCL, $L=4$','SCL, $L=8$','SCL, $L=16$','PCM-SCL, $L=2$','PCM-SCL, $L=4$','PCM-SCL, $L=8$','interpreter','Latex');
grid on


%%
clear
clc
load('Reduction_rate.mat');
plot(1:.5:4,Reduction_rate,'b-*');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('Reduction Rate','interpreter','Latex');
axis([1 4 0.45 0.7]);
grid on


%%
clear
clc
load('Additional_decoding_rate.mat')
SNR = [1 2 3 4 4.5];
semilogy(SNR,P_B,'k--',SNR,Additional_rate,'k-o',SNR,Additional_success_rate,'k-*');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('Rate','interpreter','Latex')
legend('$\mathrm{P}_\mathrm{B}$','额外解码率','额外解码成功率','interpreter','Latex');
axis([2 4.5 1e-03 1]);
grid on

%%
clear
clc
load('AverageLatency.mat')
SNR = 1:.5:4;
plot(SNR,latency_IS,'k-*',SNR,latency_LLI,'r-o')
xlabel('$E_b/N_0$ (dB)','interpreter','Latex')
ylabel('Average Latency (clock cycles)','interpreter','Latex')
axis([1 4 150 750]);
legend('IS PCM-SC-2 decoder','LLI PCM-SC-2 decoder');
grid on


%%
clear
clc
load('Block3_PCM_K140_CRC12_Kp24_R4688.mat');
figure
semilogy(SNR_SC,perSC,'k-+',SNR_SC,perBP,'b-*',SNR_SCL,perSCL2,'m-o',SNR_SCL,perSCL4,'c-^');
hold on
semilogy(SNR_PCM,perPCM_SC_3,'r-p');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('BER','interpreter','Latex');
axis([1 4.5 1.0e-4 1])
legend('SC','BP','SCL, $L=2$','SCL, $L=4$','PCM-SC-3','interpreter','Latex');
grid on

%%
clear
clc
load('bound.mat')

SNR = [1 2 3 4 4.5];

alpha_b = 0.38;
alpha_c = 6.9;
perNewMin = (1+alpha_b).*P_B.^2 - alpha_b.*P_B.^3;
perNewMax = (1+alpha_c).*P_B.^2 - alpha_c.*P_B.^3;

figure
semilogy(SNR,perPCM_SC,'r-+',SNR,perNewMax,'k-*',SNR,perNewMin,'k-o')
xlabel('$E_b/N_0$ (dB)','interpreter','Latex')
ylabel('PER','interpreter','Latex')
legend('PCM-SC-2','$P_{new}^{upper}$','$P_{new}^{lower}$','interpreter','Latex')
axis([1 4.5 1.0e-04 1.4])
grid on
