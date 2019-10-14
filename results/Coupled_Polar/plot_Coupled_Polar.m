%%%%%%%%%%%%%%%%%%%%%%
%%
clear
clc

load('Coupled_Polar_N2048_R5_Rc15.mat')

figure
semilogy(SNR, perCoupled, 'k-*');
xlabel('$E_b/N_0$','interpreter','Latex')
ylabel('PER','interpreter','Latex')
axis([1 2.5 1.0e-05 1])
legend('Coupled Polar, $L=32$', 'interpreter', 'Latex')
grid on