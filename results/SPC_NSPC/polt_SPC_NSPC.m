%%%%%%%%%%%%%%%%%%%%%%
%%
load('SPC_NSPC.mat')
figure
semilogy(SNR, berSPC, 'k-*', SNR, berNSPC, 'r-+');
xlabel('$E_b/N_0$','interpreter','Latex')
ylabel('BER','interpreter','Latex')
% axis()
legend('SPC','NSPC', 'interpreter', 'Latex')
grid on

figure
semilogy(SNR, perSPC, 'k-*', SNR, perNSPC, 'r-+');
xlabel('$E_b/N_0$','interpreter','Latex')
ylabel('PER','interpreter','Latex')
% axis()
legend('SPC','NSPC', 'interpreter', 'Latex')
grid on