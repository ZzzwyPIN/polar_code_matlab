clc
clear
addpath('../GA/');
% 基本参数设置
n = 8;  % 比特位数
R = 0.4531;    % 码率
SNR = 4;
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
Kpure = 116;
Kinfo = 128;
K = 140;  % information bit length
k_f = N - K;

% get generate matrix
G = encoding_matrix(n);
rng('shuffle')
for i = 1:length(SNR)
    sigma = (2*esn0(i))^(-0.5);
    
    % get information bits and concatenated bits
    P = GA(sigma,N);   % load the channel information
    [~, I] = sort(P,'descend');
%     info_index = I(1:K);  % 挑选质量好的信道传输信息位
    pure_index = I(1:Kpure);
    crc_index = I(Kinfo+1:K);
    info_index = [pure_index crc_index];
    frozen_index = I(K+1:end);   % 传输冻结位的信道

    % set PER and BER counter
    PerNum = 0;
    BerNum = 0;
    iter = 0;
    while true
        iter = iter + 1;
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow perNum: %2d\tNow berNum: %2d', iter, SNR(i),PerNum,BerNum);
        source_bit = randi([0 1],1,Kinfo);
        
        
        u = zeros(1,N);
        u(info_index) = source_bit;
        
        encode_temp = rem(u*G,2);
%         encode_temp = rem(source_bit*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp = 1 - 2 * encode_temp;
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        lr_x = 2*receive_sample/(sigma^2);
        % decoding follow
        lr_u = zeros(1,N); % save send sample LR in each iteration
        lr_u(reverse_index(n,frozen_index)) = init_max;
        
        receive_bits = polarBP_decoder(n,lr_u,lr_x,max_iter,pure_index);
        
        % calculate BER and PER
        count = sum(receive_bits ~= source_bit(1:Kpure));
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end
        
        if (PerNum>=100)
            break;
        end
    end
    iterNum(i) = iter;
    perBP(i) = PerNum/iter;
    berBP(i) = BerNum/(Kpure*iter);
end