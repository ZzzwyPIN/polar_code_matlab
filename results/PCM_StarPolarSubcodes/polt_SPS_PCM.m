%%%%%%%%%%%%%%%%%%%%%%
%%
clc
clear
load('PCM_SPS_Comparison.mat')
figure
semilogy(SNR, berPCM, 'k-*', SNR, berSPS, 'r-+');
xlabel('$E_b/N_0$','interpreter','Latex')
ylabel('BER','interpreter','Latex')
% axis()
legend('PCM-SCL, $L=32$','Star Polar Subcodes', 'interpreter', 'Latex')
grid on

figure
semilogy(SNR, perPCM, 'k-*', SNR, perSPS, 'r-+');
xlabel('$E_b/N_0$','interpreter','Latex')
ylabel('PER','interpreter','Latex')
% axis()
legend('PCM-SCL, $L=32$','Star Polar Subcodes', 'interpreter', 'Latex')
grid on