clear
clc

load('N1024_R5_L8_Comparison.mat')

figure
semilogy(SNR(1:end-1),per, 'k-*',SNR(1:end-1),perCP,'r-+',SNR(1:end-1),perImp,'b-<')
hold on
semilogy(SNR,perSPC,'k-+',SNR,perImpSPC,'r-h',SNR,perImpPuncSPC,'b-s',SNR,perInter,'m-^')
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('PER','interpreter','Latex')
legend('PCM','Coupled Polar Codes','PCM (Imp)','PCM (SPC)','PCM (Imp&SPC)','PCM (Imp&SPC&Punc)','interFrame')
axis([0.5 2.5 1.0e-05 1])
grid on

figure
semilogy(SNR(1:end-1),ber, 'k-*',SNR(1:end-1),berImp,'b-<')
hold on
semilogy(SNR,berSPC,'k-+',SNR,berImpSPC,'r-h',SNR,berImpPuncSPC,'k-s',SNR,berInter,'m-^')
xlabel('$E_b/N_0$ (dB)', 'interpreter', 'Latex')
ylabel('BER','interpreter','Latex')
legend('PCM','PCM (Imp)','PCM (SPC)','PCM (Imp&SPC)','PCM (Imp&SPC&Punc)','interFrame')
axis([0.5 2.5 1.0e-06 1])
grid on