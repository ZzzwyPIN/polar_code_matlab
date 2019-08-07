%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file shows the figures for optimal PCM under the code length 1024
clear
clc

%% R = 0.5
clear;
clc;
load('PCM_N1024_Kp74_snr1t3.5_R5.mat')

figure
semilogy(SNR,perSC,'k-*',SNR,per,'b-+',SNR,perSCL2,'r-o',SNR,perSCL4,'m-<');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('PER','interpreter','Latex');
legend('SC','PCM-SC-2','SCL, $L=2$','SCL, $L=4$','interpreter','Latex')
axis([1 3.5 1e-04 1])
grid on

figure
semilogy(SNR,berSC,'k-*',SNR,ber,'b-+',SNR,berSCL2,'r-o',SNR,berSCL4,'m-<');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('BER','interpreter','Latex');
legend('SC','PCM-SC-2','SCL, $L=2$','SCL, $L=4$','interpreter','Latex')
axis([1 3.5 1e-05 1])
grid on