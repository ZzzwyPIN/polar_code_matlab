function [c_index,f_index] = coding_index(snr, coding_level, code_rate)
% This function calculates the coding index for a given transition
% probability of a BSC channel
% snr : SNR of the BSC in dB
% coding_level: coding level
% code_rate: code rate

% block length
N = 2^coding_level;
% # of bits in one block
N_i = N * code_rate;

snr_l = 10^(snr/10);
% the size of Z 
Z = zeros(1,N);
% channel transition probability
p = qfunc((snr_l).^0.5);
% initial Z
Z(1) = 2*(p*(1-p))^0.5;
Z = calculate_z(Z);
[Z_s, I] = sort(Z);
c_index = sort( I(1:N_i)); 
f_index = sort(I(N_i+1:end));
end

