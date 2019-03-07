clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.75;    % 码率
SNR = 3;

% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
N = 2^n;
K = floor(N*R);  % information bit length
k_f = N - K;

% get information bits and concatenated bits
load('Pe_N256_snr3.2_R5.mat');   % load the channel information
[Ptmp, I] = sort(P);
info_index = sort(I(1:K));  % 挑选质量好的信道传输信息位
frozen_index = sort(I(K+1:end));   % 传输冻结位的信道

% get generate matrix
G = encoding_matrix(n);
Gi = G(info_index,:);
Gf = G(frozen_index,:);
frozen_bits = zeros(1,k_f);
rng('shuffle')
for i = 1:length(SNR)
    sigma = (2*esn0(i))^(-0.5);
    % set PER and BER counter
    PerNum = 0;
    BerNum = 0;
    iter = 0;
    while (true)
        iter = iter + 1;
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow perNum: %2d\tNow berNum: %2d', iter, SNR(i),PerNum,BerNum);
        source_bit = randi([0 1],1,K);
        encode_temp = rem(source_bit*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp = (-1).^(encode_temp + 1);
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        receive_bits = polarSC_decoder(n,receive_sample,sigma,frozen_index,frozen_bits,info_index);
        
        % calculate BER and PER
        count = sum(receive_bits ~= source_bit);
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end 
        if (iter >= 1000)
            break;
        end
    end
    iterNum(i) = iter;
    perSC(i) = PerNum/iter;
    berSC(i) = BerNum/(K*iter);
end