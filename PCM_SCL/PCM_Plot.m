clc
clear

%%
load('PCM_N256_R5.mat');

figure
semilogy(SNR,berSCL8,'k-+',SNR2,ber2,'k-*',SNR4,ber4,'b-^',SNR8,ber8,'r-s',SNR16,ber16,'c-o');
xlabel('Eb/N0')
ylabel('Bit Error Rate')
axis([0 3 9.40e-07 1])
legend('SCL8','PCM\_SCL2','PCM\_SCL4','PCM\_SCL8','PCM\_SCL16')
grid on

figure
semilogy(SNR,perSCL8,'k-+',SNR2,per2,'k-*',SNR4,per4,'b-^',SNR8,per8,'r-s',SNR16,per16,'c-o');
xlabel('Eb/N0')
ylabel('Packet Error Rate')
axis([0 3 4.50e-06 1])
legend('SCL8','PCM\_SCL2','PCM\_SCL4','PCM\_SCL8','PCM\_SCL16')
grid on