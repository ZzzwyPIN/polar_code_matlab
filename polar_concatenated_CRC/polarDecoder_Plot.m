% clc;
% clear;
% %%
%R = 0.25
% load('result/final/poalr_R0.25_total_1204.mat');
%ber
SNR1 = [0 1 2 3 3.5];
figure
semilogy(SNR1,berSC,'b-*',SNR1,berBP,'k-+',SNR1,berSC_ReSC_w,'r-d',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(SNR1,berSC_ReSC_b,'r-s',SNR2,berBP_ReBP,'b-p');
xlabel('Eb/N0');
ylabel('Bit Error rate');
axis([0 4 5.0e-06 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-BP')
%per
figure
semilogy(SNR1,perSC,'b-*',SNR1,perBP,'k-+',SNR2,perSC_ReSC_w,'r-d',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
semilogy(SNR1,perSC_ReSC_b,'r-s',SNR2,perBP_ReBP,'b-p',SNR2,perSC_th,'k-h');%SNR1,perSC_,'m-p',SNR1,perSC_th,'k-h','perSC^2','perSC_{bound}'
xlabel('Eb/N0');
ylabel('Frame Error rate');
axis([0 3.5 1.0e-04 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-BP','perSC_{bound}')

%%
% R = 0.5
% clear;
% load('result/final/poalr_R0.5_total_1204.mat')
%ber
SNR = [0 1 2 3 3.5 4];
SNR1 = [0 1 2 3 4];
SNR2 = [0 1 2 3 3.2];
figure
semilogy(SNR1,berSC,'b-*',SNR1,berBP,'k-+',SNR2,berSC_ReSC_w,'r-d',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(SNR2,berSC_n3,'r-s',SNR2,berBP_ReBP,'b-p');
xlabel('Eb/N0');
ylabel('Bit Error rate');
axis([0 4 5.0e-06 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SCL8','SCL32','SC-SC_{n3}','BP-BP')
%per
figure
semilogy(SNR1,perSC,'b-*',SNR1,perBP,'k-+',SNR2,perSC_ReSC_w,'r-d',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
semilogy(SNR2,perSC_n3,'r-s',SNR2,perBP_ReBP,'b-p',SNR2,perSC_th,'k-h');%SNR1,perSC_,'m-p',SNR1,perSC_th,'k-h','perSC^2','perSC_{bound}'
xlabel('Eb/N0');
ylabel('Frame Error rate');
axis([0 3.5 1.0e-04 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SCL8','SCL32','SC-SC_{n3}','BP-BP','perSC_{bound}')

%%
% R = 0.75
% clear;
% load('result/final/poalr_R0.75_total_1204.mat')
SNR = [0 1 2 3 4];
SNR1 = [1 2 3 3.5];
%ber
figure
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR1,berSC_SC_w,'r-d',SNR1,berSCL8,'c-o',SNR1,berSCL32,'b-^',SNR1,berBP_BP,'b-p');
xlabel('Eb/N0');
ylabel('Bit Error rate');
% axis([0 4 5.0e-06 1])
grid on
legend('SC','BP','SC-SC','SCL8','SCL32','BP-BP') 
%per
figure
semilogy(SNR,perSC,'b-*',SNR,perBP,'k-+',SNR1,perSC_SC_w,'r-d',SNR1,perSCL8,'c-o',SNR1,perSCL32,'b-^',SNR1,perBP_BP,'b-p');
xlabel('Eb/N0');
ylabel('Bit Error rate');
% axis([0 4 5.0e-06 1])
grid on
legend('SC','BP','SC-SC','SCL8','SCL32','BP-BP')