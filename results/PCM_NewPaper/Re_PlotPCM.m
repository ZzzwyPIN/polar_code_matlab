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


%%
clear
clc

load('PCM_SC_BP.mat')
figure
h1 = semilogy(SNR_,perSC,'-*','LineWidth',1.2);
h1.Color = [0.8549    0.6471    0.1255];
hold on
h2 = semilogy(SNR_,perBP,'-s','LineWidth',1.2);
h2.Color = 'k';
h3 = semilogy(SNR,PCM_SC_2,'-o','LineWidth',1.2);
h3.Color = 'r';
h4 = semilogy(SNR,PCM_BP_2,'-^','LineWidth',1.2);
h4.Color = [1 0 1];
h5 = semilogy(SNR,perSCL2,'-h','LineWidth',1.2);
h5.Color = [0.4 0.4 0.4];
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('PER','interpreter','Latex')
legend('SC','BP','PCM-SC-$2$','PCM-BP-$2$','SCL, $L=2$','interpreter','Latex')
axis([2 3.8 9.5339e-06 0.2])
% axis tight
grid on


figure
h1 = semilogy(SNR_,perSC,'-*','LineWidth',1.2);
h1.Color = [0.8549    0.6471    0.1255];
hold on
h2 = semilogy(SNR_,perBP,'-s','LineWidth',1.2);
h2.Color = 'k';
h3 = semilogy(SNR,PCM_SC_3,'-o','LineWidth',1.2);
h3.Color = 'r';
h5 = semilogy(SNR,perSCL2,'-h','LineWidth',1.2);
h5.Color = [0.4 0.4 0.4];
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('PER','interpreter','Latex')
legend('SC','BP','PCM-SC-$3$','SCL, $L=2$','interpreter','Latex')
% axis([2 3.8 3.0e-05 0.2])
axis tight
grid on


%%
clear
clc
load('Plot_AdditionalRate.mat');

figure
semilogy(SNR,PB,'k-*',SNR,Pb,'r-+','LineWidth',1.2);
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('PER','interpreter','Latex')
legend('PB','Additional successful decoding rate')
grid on

figure
plot(SNR,Rate,'k-*','LineWidth',1.2);
xlabel('$E_b/N_0$ (dB)', 'interpreter','Latex')
ylabel('Ratio','interpreter','Latex')
axis([2 3.6 0 10])
grid on


%%
clear
clc

load('PCM_SCL.mat')
figure
h1 = semilogy(SNR,PCM_SCL2,'-*','LineWidth',1.2);
h1.Color = [0.8549    0.6471    0.1255];
hold on
h2 = semilogy(SNR,PCM_SCL4,'-s','LineWidth',1.2);
h2.Color = [1 0 1];
h3 = semilogy(SNR(1:end-1),PCM_SCL8,'-o','LineWidth',1.2);
h3.Color = 'r';
h4 = semilogy(SNR,SCL4,'-^','LineWidth',1.2);
h4.Color = 'k';
h5 = semilogy(SNR,SCL8,'-p','LineWidth',1.2);
h5.Color = [0.4 0.4 0.4];
h6 = semilogy(SNR,SCL16,'-h','LineWidth',1.2);
h6.Color = [0.5451    0.3961    0.0314];
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('PER','interpreter','Latex')
legend('PCM\_SCL, $L=2$','PCM\_SCL, $L=4$','PCM\_SCL, $L=8$','SCL, $L=4$','SCL, $L=8$','SCL, $L=16$','interpreter','Latex')
axis([1 3 6.6613e-07 1])
% axis tight
grid on

%%
clear
clc

load('DerivedPER.mat')

figure
semilogy(SNR,Pnew,'k-*',SNR,Pnew_,'k-+',SNR,PCM_SC_2,'r-o','LineWidth',1.2)
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('PER','interpreter','Latex')
legend('Upper Bound','Lower Bound','PCM-SC-2','interpreter','Latex')
axis tight
grid on


%%
clear
clc

load('Additional_decoding.mat')
semilogy(SNR,PB,'k--',SNR,Add_Rate,'r-*',SNR,Add_Success_Rate,'k-o','LineWidth',1.2)
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('Rate','interpreter','Latex')
legend('PB','Additional Decoding Rate','Additional Success Rate','interpreter','Latex')
axis tight
grid on

%%
clear
clc
load('N1024_R5_L8_Comparison.mat')

figure
h1 = semilogy(SNR,perInter,'-+','LineWidth',1.2);
h1.Color = [0.4 0.4 0.4];
hold on
h2 = semilogy(SNR(1:end-1),perCP,'-*','LineWidth',1.2);
h2.Color = 'k';
h3 = semilogy(SNR,per,'-o','LineWidth',1.2);
h3.Color = 'r';
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('PER','interpreter','Latex')
legend('[12]','[13]','PCM')
axis([0.5 2.5 1.0e-05 1])
grid on