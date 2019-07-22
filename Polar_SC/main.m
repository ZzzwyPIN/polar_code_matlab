clc
clear
addpath('../GA/');
% 基本参数设置
n = 8;  % 比特位数
R = 0.4531;    % 码率
SNR = [2 3 4 4.5 5];

% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
N = 2^n;
Kpure = 116;
Kinfo = 128;
K = 140;  % information bit length
k_f = N - K;


% get generate matrix
G = encoding_matrix(n);
frozen_bits = zeros(1,k_f);

rng('shuffle')
for i = 1:length(SNR)
    sigma = (2*esn0(i))^(-0.5);
    
    % get information bits and concatenated bits
%     load('Pe_N256_snr3.2_R5.mat');   % load the channel information
    P = GA(sigma,N);
    [~, I] = sort(P,'descend');
    pure_index = I(1:Kpure);
    crc_index = I(Kinfo+1:K);
    info_index = [pure_index crc_index];
    frozen_index = I(K+1:end);   % 传输冻结位的信道

   
    % set PER and BER counter
    PerNum = 0;
    BerNum = 0;
    iter = 0;
    while (true)
        iter = iter + 1;
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow perNum: %2d\tNow berNum: %2d', iter, SNR(i),PerNum,BerNum);
        source_bit = randi([0 1],1,Kinfo);
        
        u = zeros(1,N);
        u(info_index) = source_bit;
        encode_temp = rem(u*G,2);
    
        % bpsk modulation
        encode_temp = 1-2*encode_temp;
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        receive_bits = polarSC_decoder(n,receive_sample,sigma,frozen_index,frozen_bits,pure_index);
        
        % calculate BER and PER
        count = sum(receive_bits ~= source_bit(1:Kpure));
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end 
        if ( PerNum >=100)
            break;
        end
    end
    iterNum(i) = iter;
    perSC(i) = PerNum/iter;
    berSC(i) = BerNum/(Kpure*iter);
end
