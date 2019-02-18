function SNR = SNR_calculate(R, SNRindB)
esn0 = 10^(SNRindB/10);
snr = esn0/R;
SNR = 10*log10(snr);