clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.2188;    % 码率
SNR = -3:2;
block_num = 10000;

% 参数计算
snr = 10.^(SNR/10);
N = 2^n;
K = floor(N*R);  % information bit length
k_f = N - K;

% get information bits and concatenated bits
load('Pe_snr3p0db_2048_n_8.mat');   % load the channel information
[Ptmp, I] = sort(P);
Info_index = sort(I(K:-1:1));  % 挑选质量好的信道传输信息位
Frozen_index = sort(I(end:-1:K+1));   % 传输冻结位的信道

% get generate matrix
G = encoding_matrix(n);
Gi = G(Info_index,:);
Gf = G(Frozen_index,:);
frozen_bits = zeros(1,k_f);
rng('shuffle')
for i = 1:length(SNR)
    sigma = (1/snr(i))^0.5;
    % set PER and BER counter
    PerNum = 0;
    BerNum = 0;
    for iter = 1:block_num
        fprintf('\nNow iter: %2d\tNow SNR: %d', iter, SNR(i));
        source_bit = randi([0 1],1,K);
        encode_temp = rem(source_bit*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp = (-1).^(encode_temp + 1);
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        receive_bits = polarSC_decoder(n,receive_sample,snr(i),Frozen_index,frozen_bits,Info_index);
        
        % calculate BER and PER
        count = sum(receive_bits ~= source_bit);
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end 
    end
    perSC(i) = PerNum/block_num;
    berSC(i) = BerNum/(K*block_num);
end