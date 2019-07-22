clc
clear
%% The PER of PCM for variable Kp under the Eb/N0=4 dB. N=256
clear
clc
load('Optimal_PCM_N256_Kp2t2t64_snr4.mat');
figure
semilogy(Kp,per,'k-*');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('PER','interpreter','Latex');
axis([0 64 1e-04 1e-02]);
grid on;

figure
plot(Kp,Rate,'r-o');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('Rate','interpreter','Latex');
% axis([0 64 3e-04 1e-02]);
grid on;


%% The PER of PCM for variable Kp under the Eb/N0=3.5 dB. N=512
clear
clc
load('Optimal_PCM_N512_Kp2t2t128_snr3d5.mat');
figure
semilogy(Kp,per,'k-*');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('PER','interpreter','Latex');
% axis([0 64 3e-04 1e-02]);
grid on;

figure
plot(Kp,Rate,'r-o');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('Rate','interpreter','Latex');
% axis([0 64 3e-04 1e-02]);
grid on;

%% The PER of PCM for variable Kp under the Eb/N0=3 dB. N=1024
clear
clc
load('Optimal_PCM_N1024_Kp2t4t256_snr3.mat');
figure
semilogy(Kp,per,'k-*');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('PER','interpreter','Latex');
% axis([0 64 3e-04 1e-02]);
grid on;

figure
plot(Kp,Rate,'r-o');
xlabel('$K_{\mathrm{p}}$','interpreter','Latex');
ylabel('Rate','interpreter','Latex');
% axis([0 64 3e-04 1e-02]);
grid on;








%% plot the PER of variable Kp
clc
clear
SNR = [1 2 3 4 4.5];
load('Optimal_Kp.mat')

semilogy(SNR,perSC_Kp10,'k-+',SNR,perSC_Kp24,'r-s',SNR,perSC_Kp32,'c-^');
hold on
semilogy(SNR,perPCM_Kp10,'b-*',SNR,perPCM_Kp24,'m-h',SNR,perPCM_Kp32,'b-o');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('PER','interpreter','Latex');
legend('SC, $K_p=10$','SC, $K_p=24$','SC, $K_p=32$','PCM, $K_p=10$','PCM, $K_p=24$','PCM, $K_p=32$','interpreter','Latex');
% axis([])
grid on;

%% plot the rate of perSC and perPCM for different Kp
clc
clear
SNR = [1 2 3 4 4.5];
Kp = 2:2:32;
load('PER_Rate_Kp2to32_SNR1to4dot5.mat')
figure
plot(SNR,rate(5,:),'c-h');
hold on
plot(SNR,rate(6,:),'k-+',SNR,rate(7,:),'b-*',SNR,rate(8,:),'r-o',SNR,rate(9,:),'m-s',SNR,rate(10,:),'c-v');
plot(SNR,rate(11,:),'k-^',SNR,rate(12,:),'b-<',SNR,rate(13,:),'r->',SNR,rate(14,:),'m--s',SNR,rate(15,:),'c--h');
plot(SNR,rate(16,:),'k--^')
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('Rate','interpreter','Latex');
legend('$K_p=10$','$K_p=12$','$K_p=14$','$K_p=16$','$K_p=18$','$K_p=20$','$K_p=22$','$K_p=24$','$K_p=26$','$K_p=28$','$K_p=30$','$K_p=32$','interpreter','Latex');
axis([1 4.5 0 10])
% axis([1 3 0.5 2])
grid on


figure
plot(Kp,rate(:,1),'k-+',Kp,rate(:,2),'b-*',Kp,rate(:,3),'r-o',Kp,rate(:,4),'m-s',Kp,rate(:,5),'c-h');
xlabel('$E_b/N_0$ (dB)','interpreter','Latex');
ylabel('Rate','interpreter','Latex');
legend('$SNR=1$ (dB)','$SNR=2$ (dB)','$SNR=3$ (dB)','$SNR=4$ (dB)','$SNR=4.5$ (dB)','interpreter','Latex');
axis([0 32 0 10])
grid on












