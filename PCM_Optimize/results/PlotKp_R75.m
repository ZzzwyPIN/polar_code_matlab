clc
clear
%% The PER of PCM for variable Kp under the Eb/N0=5 dB. N=256 R=0.75
clear
clc
load('Optimal_PCM_N256_Kp2t4t96_snr5_R75.mat');
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

%% The PER of PCM for variable Kp under the Eb/N0=4.5 dB. N=512 R=0.75
clear
clc
load('Optimal_PCM_N512_Kp2t4t192_snr4d5_R75.mat');
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


%% The PER of PCM for variable Kp under the Eb/N0=4 dB. N=1024 R=0.75
clear
clc
load('Optimal_PCM_N1024_Kp4t4t384_snr4_R75.mat');
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