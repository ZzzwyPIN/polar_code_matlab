clc;
clear;
%%
%R = 0.25
load('result/final/poalr_R0.25_total_1204.mat');
esn0 = [-0.4375 0 0.4375 0.8750 1.3125 1.7500 2.1875];
esn0_ = [0 0.4375 0.8750 1.3125 1.5313 1.7500];
%ber
figure
semilogy(esn0,berSC,'b-*',esn0,berBP,'k-+',esn0,berSC_ReSC_w,'r-d',esn0,berSC_ReBP,'m-p',esn0,berSCL8,'c-o',esn0,berSCL32,'b-^');
hold on
semilogy(0:4,berSC_ReSC_b,'r-s',0:4,berBP_ReSC,'k-h',0:4,berBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([0 4 6.1e-06 0.4])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')

%per
figure
semilogy(esn0,perSC,'b-*',esn0,perBP,'k-+',esn0,perSC_ReSC_w,'r-d',esn0,perSC_ReBP,'m-p',esn0,perSCL8,'c-o',esn0,perSCL32,'b-^');
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
esn0 = [-0.4375 0 0.4375 0.8750 1.3125 1.7500 2.1875];
esn0_ = [0 0.4375 0.8750 1.3125 1.5313 1.7500];
figure
semilogy(esn0_,berSC,'b-*',esn0_,berBP,'k-+',esn0_,berSC_ReSC_w,'r-d',esn0_,berSC_ReBP,'m-p',esn0,berSCL8,'c-o',esn0,berSCL32,'b-^');
hold on
semilogy(esn0_,berSC_ReSC_b,'r-s',esn0_,berBP_ReSC,'k-h',esn0_,berBP_ReBP,'b-p');
xlabel('Eb/N0');
ylabel('Bit Error rate');
axis([0 1.75 5.0e-06 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')
%per
figure
semilogy(esn0_,perSC,'b-*',esn0_,perBP,'k-+',esn0_,perSC_ReSC_w,'r-d',esn0_,perSC_ReBP,'m-p',esn0,perSCL8,'c-o',esn0,perSCL32,'b-^');
hold on
semilogy(esn0_,perSC_ReSC_b,'r-s',esn0_,perBP_ReSC,'k-h',esn0_,perBP_ReBP,'b-p');
xlabel('Eb/N0');
ylabel('Frame Error rate');
axis([0 1.75 2.8e-05 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')

%%
%R = 0.75
clear;
load('result/final/poalr_R0.75_total_1204.mat')
esn0 = [-0.4375 0 0.4375 0.8750 1.3125 1.7500 2.1875];
esn0_ = [0 0.4375 0.8750 1.3125 1.5313 1.7500];
%ber
figure
semilogy(esn0,berSC,'b-*',esn0,berBP,'k-+',esn0,berSC_ReSC_w,'r-d',esn0,berSC_ReBP,'m-p',esn0,berSCL8,'c-o',esn0,berSCL32,'b-^');
hold on
semilogy(0:4,berSC_ReSC_b,'r-s',0:4,berBP_ReSC,'k-h',0:4,berBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([1 4 1.0e-05 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')
%per
figure
semilogy(esn0,perSC,'b-*',esn0,perBP,'k-+',esn0,perSC_ReSC_w,'r-d',esn0,perSC_ReBP,'m-p',esn0,perSCL8,'c-o',esn0,perSCL32,'b-^');
hold on
semilogy(0:4,perSC_ReSC_b,'r-s',0:4,perBP_ReSC,'k-h',0:4,perBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Frame Error rate');
axis([1 5 6.00e-05 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')