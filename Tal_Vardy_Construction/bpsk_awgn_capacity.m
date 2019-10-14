clc; clear;
close all
x = -100:1:100;

d = 0.01;
%snr = 0.2:d:5; % linear SNR
snr = -3:0.1:15; % dB
%snr = 10^(0.1/10) 
I = zeros(1,length(snr));

for i=1:length(snr)
    % fading loss
    fs = 0;
    % amplitude
    snr_linear = 10^((snr(i)-fs)/20);
    %a = snr(i);
    a = snr_linear;
    fa = 0;
    fma = 0;
    for j=1:length(x)
        fa = fa+1/(2*pi)^0.5*exp(-(x(j)-a)^2/2)*log2(2/(1+exp(-2*a*x(j))));
        fma = fma+1/(2*pi)^0.5*exp(-(x(j)+a)^2/2)*log2(2/(1+exp(2*a*x(j))));
    end 
    I(i) = 0.5*fa + 0.5*fma;
end
%snr_db = 20*log10(snr);
snr_db = snr;
plot(snr_db,I)
xlabel('SNR (dB)')
ylabel('Capacity');

