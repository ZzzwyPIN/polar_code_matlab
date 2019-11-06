%%
clear;
clc;

load('Re_PCM.mat')

figure
h1 = semilogy(SNR,perSC(1:end-1),'-*','LineWidth',1.2);
h1.Color = [0.8549    0.6471    0.1255];
hold on
h2 = semilogy(SNR,perBP(1:end-1),'-s','LineWidth',1.2);
h2.Color = 'k';
h3 = semilogy(SNR,perSC_PCM,'-o','LineWidth',1.2);
h3.Color = 'r';
h4 = semilogy(SNR,perBP_PCM,'-^','LineWidth',1.2);
h4.Color = [1 0 1];
h5 = semilogy(SNR,perSCL2,'-h','LineWidth',1.2);
h5.Color = [0.4 0.4 0.4];
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('PER','interpreter','Latex')
legend('SC','BP','PCM-SC-$2$','PCM-BP-$2$','SCL, $L=2$','interpreter','Latex')
axis([1 3.5 1.0e-04 1])
grid on


figure
p1 = semilogy(SNR,berSC(1:end-1),'-*','LineWidth',1.2);
p1.Color = [0.8549    0.6471    0.1255];
hold on
p2 = semilogy(SNR,berBP(1:end-1),'-s','LineWidth',1.2);
p2.Color = 'k';
p3 = semilogy(SNR,berSC_PCM,'-o','LineWidth',1.2);
p3.Color = 'r';
p4 = semilogy(SNR,berBP_PCM,'-^','LineWidth',1.2);
p4.Color = [1 0 1];
p5 = semilogy(SNR,berSCL2,'-h','LineWidth',1.2);
p5.Color = [0.4 0.4 0.4];
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('BER','interpreter','Latex')
legend('SC','BP','PCM-SC-$2$','PCM-BP-$2$','SCL, $L=2$','interpreter','Latex')
axis([1 3.5 1.0e-05 1])
grid on