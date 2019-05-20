clc
clear

% 基本参数设置
n = 10;  % 比特位数
R = 0.4375;    % 码率
SNR = [1 2 3 3.5 4];
init_lr_max = 3;    % limit the max LR of the channel to be with [-3 3]
max_iter = 40;

% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
N = 2^n;
K = floor(N*R);  % information bit length
k_f = N - K;

[M_up, M_down] = index_Matrix(N);
lambda_offset = 2.^(0 : n);
llr_layer_vec = get_llr_layer(N);

% get information bits and concatenated bits
load('GA_N1024_R5_snr3.2.mat');   % load the channel information
[Ptmp, I] = sort(P,'descend');
info_index = sort(I(K:-1:1));  % 挑选质量好的信道传输信息位

% get generate matrix
frozen_bits = ones(N , 1);
frozen_bits(info_index) = 0;
rng('shuffle')
for i = 1:length(SNR)
    sigma = (2*esn0(i))^(-0.5);
    % set PER and BER counter
    PerNum = 0;
    BerNum = 0;
    iter = 0;
    while true
        iter = iter + 1;
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow perNum: %2d\tNow berNum: %2d', iter, SNR(i),PerNum,BerNum);
        u = zeros(N, 1);
        source_bit = randi([0 1],1,K);
        u(info_index) = source_bit;
        encode_temp = polar_encoder(u, lambda_offset, llr_layer_vec);
    
        % bpsk modulation
        encode_temp = 1 - 2*encode_temp;
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        llr = 2*receive_sample/(sigma^2);
        
        [receive_bits, ~, ~, ~] = BP_Decoder_LLR(info_index, frozen_bits, llr, max_iter, M_up, M_down);
        
        % calculate BER and PER
        count = sum(receive_bits' ~= source_bit);
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end
        
        if (iter >= 10000 && PerNum >= 100)
            break;
        end
    end
    iterNum(i) = iter;
    perBP(i) = PerNum/iter;
    berBP(i) = BerNum/(K*iter);
end