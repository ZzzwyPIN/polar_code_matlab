
%%
clear
clc
load('N1024_R5_L8_Comparison.mat')

figure
h1 = semilogy(SNR,perInter,'-+','LineWidth',1.2);
h1.Color = 'k';
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

% figure
% semilogy(SNR,ber,'k-*',SNR(1:end-1),berImp,'m-o','LineWidth',1.2)
% hold on
% semilogy(SNR,berSPC_Punc,'b-s',SNR,berInter,'k-^',SNR,berPunc,'r-h','LineWidth',1.2)
% xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
% ylabel('BER','interpreter','Latex')
% legend('PCM','PCM (Imp)','PCM (Imp&SPC&Punc)','interFrame','PCM (Punc)')
% axis([0.5 2.5 1.0e-06 1])
% grid on