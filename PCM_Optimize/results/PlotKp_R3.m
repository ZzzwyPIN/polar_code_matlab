clc
clear
%% The PER of PCM for variable Kp under the Eb/N0=4 dB. N=256 R=1/3
clear
clc
load('Optimal_PCM_N256_Kp2t2t42_snr4_R3.mat');
figure
semilogy(Kp,per,'k-*');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('PER','interpreter','Latex');
% axis([0 64 1e-04 1e-02]);
grid on;

figure
plot(Kp,Rate,'r-o');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('Rate','interpreter','Latex');
% axis([0 64 3e-04 1e-02]);
grid on;

%% The PER of PCM for variable Kp under the Eb/N0=3.5 dB. N=512 R=1/3
clear
clc
load('Optimal_PCM_N512_Kp2t2t84_snr3d5_R3.mat');
figure
semilogy(Kp,per,'k-*');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('PER','interpreter','Latex');
% axis([0 64 1e-04 1e-02]);
grid on;

figure
plot(Kp,Rate,'r-o');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('Rate','interpreter','Latex');
% axis([0 64 3e-04 1e-02]);
grid on;


%% The PER of PCM for variable Kp under the Eb/N0=3 dB. N=1024 R=1/3
clear
clc
load('Optimal_PCM_N1024_Kp2t4t172_snr3_R3.mat');
figure
semilogy(Kp,per,'k-*');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('PER','interpreter','Latex');
% axis([0 64 1e-04 1e-02]);
grid on;

figure
plot(Kp,Rate,'r-o');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('Rate','interpreter','Latex');
% axis([0 64 3e-04 1e-02]);
grid on;