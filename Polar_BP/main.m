clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.5;    % 码率
SNR = [0 1 2 3 4];
init_lr_max = 3;    % limit the max LR of the channel to be with [-3 3]
max_iter = 40;

% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
init_max = init_lr_max * n;
if init_max > 30
    init_max = 30;
end
N = 2^n;
K = floor(N*R);  % information bit length
k_f = N - K;

% get information bits and concatenated bits
load('Pe_N256_snr3.2_R5.mat');   % load the channel information
[Ptmp, I] = sort(P);
info_index = sort(I(K:-1:1));  % 挑选质量好的信道传输信息位
frozen_index = sort(I(end:-1:K+1));   % 传输冻结位的信道

% get generate matrix
G = encoding_matrix(n);
Gi = G(info_index,:);
Gf = G(frozen_index,:);
% frozen_bits = randi([0 1],1,k_f);
frozen_bits = zeros(1,k_f);
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
        source_bit = randi([0 1],1,K);
        encode_temp = rem(source_bit*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp = (-1).^(encode_temp + 1);
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        lr_x = -2*receive_sample/(sigma^2);
        % decoding follow
        lr_u = zeros(1,N); % save send sample LR in each iteration
        lr_u(reverse_index(n,frozen_index)) = init_max;
%         frozen_index_0 = find(frozen_bits == 0);
%         frozen_index_1 = find(frozen_bits == 1);
%         lr_u(reverse_index(n,frozen_index(frozen_index_0))) = init_max;
%         lr_u(reverse_index(n,frozen_index(frozen_index_1))) = -init_max;
        
        receive_bits = polarBP_decoder(n,lr_u,lr_x,max_iter,info_index);
        
        % calculate BER and PER
        count = sum(receive_bits ~= source_bit);
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end
        
        if (iter >= 100)
            break;
        end
    end
    iterNum(i) = iter;
    perBP(i) = PerNum/iter;
    berBP(i) = BerNum/(K*iter);
end