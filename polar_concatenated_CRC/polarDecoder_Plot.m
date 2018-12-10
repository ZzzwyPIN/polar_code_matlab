clc;
clear;
%%
%R = 0.25
load('result/final/poalr_R0.25_total_1204.mat');

%ber
figure
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR,berSC_ReSC_w,'r-d',SNR,berSC_ReBP,'m-p',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(0:4,berSC_ReSC_b,'r-s',0:4,berBP_ReSC,'k-h',0:4,berBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([0 4 6.1e-06 0.4])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')

%per
figure
semilogy(SNR,perSC,'b-*',SNR,perBP,'k-+',SNR,perSC_ReSC_w,'r-d',SNR,perSC_ReBP,'m-p',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
semilogy(0:4,perSC_ReSC_b,'r-s',0:4,perBP_ReSC,'k-h',0:4,perBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Frame Error rate');
axis([0 4 2.2e-05 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')

%%
%R = 0.5
clear;
load('result/final/poalr_R0.5_total_1204.mat')
%ber
SNR = -1:5;
SNR1 = [0 1 2 3 3.2 3.4 3.6 3.8 4];
SNR2 = [0 1 2 2.5 3 3.2 3.4 3.6 3.8 4];
figure
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR1,berSC_ReSC_w,'r-d',SNR1,berSC_ReBP,'m-p',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(SNR2,berSC_ReSC_b,'r-s',SNR1,berBP_ReSC,'k-h',SNR1,berBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([0 4 9.4e-06 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')
%per
figure
semilogy(SNR,perSC,'b-*',SNR,perBP,'k-+',SNR1,perSC_ReSC_w,'r-d',SNR1,perSC_ReBP,'m-p',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
semilogy(SNR2,perSC_ReSC_b,'r-s',SNR1,perBP_ReSC,'k-h',SNR1,perBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Frame Error rate');
axis([0 4 4.4e-05 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')

%%
%R = 0.75
clear;
load('result/final/poalr_R0.75_total_1204.mat')
%ber
figure
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR,berSC_ReSC_w,'r-d',SNR,berSC_ReBP,'m-p',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(0:4,berSC_ReSC_b,'r-s',0:4,berBP_ReSC,'k-h',0:4,berBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([1 4 1.0e-05 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')
%per
figure
semilogy(SNR,perSC,'b-*',SNR,perBP,'k-+',SNR,perSC_ReSC_w,'r-d',SNR,perSC_ReBP,'m-p',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
semilogy(0:4,perSC_ReSC_b,'r-s',0:4,perBP_ReSC,'k-h',0:4,perBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Frame Error rate');
axis([1 5 6.00e-05 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')