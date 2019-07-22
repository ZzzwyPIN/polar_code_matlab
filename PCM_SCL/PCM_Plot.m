clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N = 256 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% R = 0.75
load('results/PCM_N256_R75_0524.mat');
figure
semilogy(SNR_SCL2,berSCL2,'k-+',SNR_SCL4,berSCL4,'b-d',SNR_SCL8,berSCL8,'r-p',SNR_SCL16,berSCL16,'c-h',SNR_SCL32,berSCL32,'r--+');
hold on
semilogy(SNR_PCMSCL2,berPCMSCL2,'k-*',SNR_PCMSCL4,berPCMSCL4,'b-^',SNR_PCMSCL8,berPCMSCL8,'k-s',SNR_PCMSCL16,berPCMSCL16,'r-o');
xlabel('Eb/N0')
ylabel('Bit Error Rate')
axis([1 4 2.50e-06 1])
legend('SCL2','SCL4','SCL8','SCL16','SCL32','PCM\_SCL2','PCM\_SCL4','PCM\_SCL8','PCM\_SCL16')
grid on

figure
semilogy(SNR_SCL2,perSCL2,'k-+',SNR_SCL4,perSCL4,'b-d',SNR_SCL8,perSCL8,'r-p',SNR_SCL16,perSCL16,'c-h',SNR_SCL32,perSCL32,'r--+');
hold on
semilogy(SNR_PCMSCL2,perPCMSCL2,'k-*',SNR_PCMSCL4,perPCMSCL4,'b-^',SNR_PCMSCL8,perPCMSCL8,'k-s',SNR_PCMSCL16,perPCMSCL16,'r-o');
xlabel('Eb/N0')
ylabel('Packet Error Rate')
axis([1 4 1.0e-05 1])
legend('SCL2','SCL4','SCL8','SCL16','SCL32','PCM\_SCL2','PCM\_SCL4','PCM\_SCL8','PCM\_SCL16')
grid on

%%
clear;
% R = 0.5
load('results/PCM_N256_R5_0518.mat');

figure
semilogy(SNR_SCL2,berSCL2,'k-+',SNR_SCL4,berSCL4,'b-d',SNR_SCL8,berSCL8,'r-p',SNR_SCL16,berSCL16,'c-h',SNR_SCL32,berSCL32,'r--+');
hold on
semilogy(SNR_PCMSCL2,berPCMSCL2,'k-*',SNR_PCMSCL4,berPCMSCL4,'b-^',SNR_PCMSCL8,berPCMSCL8,'r-s',SNR_PCMSCL16,berPCMSCL16,'g-o');
xlabel('Eb/N0')
ylabel('Bit Error Rate')
axis([0 3 9.40e-07 1])
legend('SCL2','SCL4','SCL8','SCL16','SCL32','PCM\_SCL2','PCM\_SCL4','PCM\_SCL8','PCM\_SCL16')
grid on

figure
semilogy(SNR_SCL2,perSCL2,'k-+',SNR_SCL4,perSCL4,'b-d',SNR_SCL8,perSCL8,'r-p',SNR_SCL16,perSCL16,'c-h',SNR_SCL32,perSCL32,'r--+');
hold on
semilogy(SNR_PCMSCL2,perPCMSCL2,'k-*',SNR_PCMSCL4,perPCMSCL4,'b-^',SNR_PCMSCL8,perPCMSCL8,'r-s',SNR_PCMSCL16,perPCMSCL16,'g-o');
xlabel('Eb/N0')
ylabel('Packet Error Rate')
axis([0 3 4.50e-06 1])
legend('SCL2','SCL4','SCL8','SCL16','SCL32','PCM\_SCL2','PCM\_SCL4','PCM\_SCL8','PCM\_SCL16')
grid on

%%
clear;
clc;
% R = 0.25
load('results/PCM_N256_R25_0522.mat');

figure
semilogy(SNR_SCL2,berSCL2,'k-+',SNR_SCL4,berSCL4,'b-d',SNR_SCL8,berSCL8,'r-p',SNR_SCL16,berSCL16,'c-h',SNR_SCL32,berSCL32,'r--+');
hold on
semilogy(SNR_PCMSCL2,berPCMSCL2,'k-*',SNR_PCMSCL4,berPCMSCL4,'b-^',SNR_PCMSCL8,berPCMSCL8,'r-s',SNR_PCMSCL16,berPCMSCL16,'g-o');
xlabel('Eb/N0')
ylabel('Bit Error Rate')
% axis([0 2 9.40e-07 1])
legend('SCL2','SCL4','SCL8','SCL16','SCL32','PCM\_SCL2','PCM\_SCL4','PCM\_SCL8','PCM\_SCL16')
grid on

figure
semilogy(SNR_SCL2,perSCL2,'k-+',SNR_SCL4,perSCL4,'b-d',SNR_SCL8,perSCL8,'r-p',SNR_SCL16,perSCL16,'c-h',SNR_SCL32,perSCL32,'r--+');
hold on
semilogy(SNR_PCMSCL2,perPCMSCL2,'k-*',SNR_PCMSCL4,perPCMSCL4,'b-^',SNR_PCMSCL8,perPCMSCL8,'r-s',SNR_PCMSCL16,perPCMSCL16,'g-o');
xlabel('Eb/N0')
ylabel('Packet Error Rate')
axis([0 3.5 6.0e-06 1])
legend('SCL2','SCL4','SCL8','SCL16','SCL32','PCM\_SCL2','PCM\_SCL4','PCM\_SCL8','PCM\_SCL16')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N=1024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%