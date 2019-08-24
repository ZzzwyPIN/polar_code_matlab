%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file shows the figures for optimal PCM under the code length 512
%% R = 0.5
clear;
clc;
load('PCM_N512_Kp42_snr1t4_R5.mat')

figure
semilogy(SNR,perSC,'k-*',SNR,perSCL2,'r-o',SNR,perSCL4,'m-<',SNR,per,'b-+',SNR,perPCM_Kp20,'c-s',SNR,perPCM_Kp64,'m-h');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('PER','interpreter','Latex');
legend('SC','SCL, $L=2$','SCL, $L=4$','PCM, $K_{\mathrm{p}}=42$','PCM, $K_{\mathrm{p}}=20$','PCM, $K_{\mathrm{p}}=64$','interpreter','Latex');
axis([1 4 1e-04 1])
grid on

figure
semilogy(SNR,berSC,'k-*',SNR,berSCL2,'r-o',SNR,berSCL4,'m-<',SNR,ber,'b-+',SNR,berPCM_Kp20,'c-s',SNR,berPCM_Kp64,'m-h');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('BER','interpreter','Latex');
legend('SC','SCL, $L=2$','SCL, $L=4$','PCM, $K_{\mathrm{p}}=42$','PCM, $K_{\mathrm{p}}=20$','PCM, $K_{\mathrm{p}}=64$','interpreter','Latex');
axis([1 4 1e-05 1])
grid on