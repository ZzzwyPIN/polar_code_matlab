function [c_index,f_index] = coding_index(snr, coding_level, code_rate, era, bec_flag)
% This function calculates the coding index for a given transition
% probability of a BSC channel
% snr : SNR of the BSC in dB
% coding_level: coding level
% code_rate: code rate
% era:
% bec_flag: The flag of BEC channel

% block length
N = 2^coding_level;
% # of bits in one block
N_i = floor(N * code_rate);
snr_l = 10^(snr/10);
% the size of Z 
Z = zeros(1,N);
% channel transition probability
p = qfunc((snr_l).^0.5);
% initial Z
if bec_flag
    Z(1) = era;
else
    Z(1) = 2*(p*(1-p))^0.5;
end
Z = calculate_z(Z);
[~, I] = sort(Z);
c_index = sort( I(1:N_i)); 
f_index = sort(I(N_i+1:end));
end

