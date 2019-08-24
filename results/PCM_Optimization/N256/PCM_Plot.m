
%%
clear
clc
load('PCM_K140_CRC12_Kp24_R4531.mat');

figure
semilogy(SNR,berSC,'b-+',SNR_SCL,berSCL2,'r-s',SNR_SCL,berSCL4,'m-^');%
hold on
semilogy(SNR_PCM(1:4),berPCM_SC,'k-o',SNR,berPCM_Kp12,'k-*',SNR,berPCM_Kp48,'b-p');%
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('BER','interpreter','Latex');
axis([1 4.5 1.0e-5 1])
legend('SC','SCL, $L=2$','SCL, $L=4$','PCM, $K_{\mathrm{p}}=24$','PCM, $K_{\mathrm{p}}=12$','PCM, $K_{\mathrm{p}}=48$','interpreter','Latex');
grid on

figure
semilogy(SNR,perSC,'b-+',SNR_SCL,perSCL2,'r-s',SNR_SCL,perSCL4,'m-^');%
hold on
semilogy(SNR_PCM,perPCM_SC,'r-o',SNR,perPCM_Kp12,'k-*',SNR,perPCM_Kp48,'m-p');%
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('PER','interpreter','Latex');
axis([1 4.5 1.0e-4 1])
legend('SC','SCL, $L=2$','SCL, $L=4$','PCM, $K_{\mathrm{p}}=24$','PCM, $K_{\mathrm{p}}=12$','PCM, $K_{\mathrm{p}}=48$','interpreter','Latex');
grid on