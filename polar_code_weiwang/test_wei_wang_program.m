% LL, Sep. 20, 2016
% test Wei Wang's program on polar codes
% info indices are selected using Tal-Vardy algorithm
clc; clear; close all;
n=8;
N=2^n;
R=0.5;

% 1000 blocks are run
SNR = 1:0.5:3; % in dB
ber = [0.1443 0.0770 0.0314 0.0089 0.0033];
per = [0.5140 0.2830 0.1250 0.0410 0.0150];

figure
semilogy(SNR, ber, 'k-o', SNR, per, 'r--');
xlabel('SNR (dB)');
legend('BER','PER');
