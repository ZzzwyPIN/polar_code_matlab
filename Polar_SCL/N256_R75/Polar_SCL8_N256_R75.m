clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.6563;    % 码率
Ng = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
L = 8;   %SCL List

SNR = [0 1 2 3 3.5];
% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
N = 2^n;

lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

K = floor(N*R);  % information bit length
Kp = floor(N*R*0.25);  % Cascaded decoding length
k_f = N-K;% frozen_bits length


filename = 'Pe_N256_snr3.2_R5.mat';
% get information bits and concatenated bits
load(filename);   % load the channel information
[Ptmp, I] = sort(P);
info_index = sort(I(1:K));  % 挑选质量好的信道传输信息位
info_without_crc = info_index(1:K-Ng);  %得到K_{info}个信息位信道
frozen_index = sort(I(K+1:end));   % 传输冻结位的信道

frozen_bits = ones(N,1);
frozen_bits(info_index) = 0;
rng('shuffle');
for i = 1:length(SNR)
    
    sigma = (2*esn0(i))^(-0.5);
    % set PER and BER counter
    PerNum = 0;
    BerNum = 0;

    iter = 0;

    
    while true 
        
        iter = iter + 1;
        % reset the frozen bits and mutual bits
        
        
        
       
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum: %2d\tNow Error Bits: %2d', iter, SNR(i),PerNum,BerNum);
        source_bit = rand(1,K-Ng)>0.5;

        source_crc_bit = crcadd(source_bit,poly);
        u = zeros(N, 1);

        u(info_index) = source_crc_bit;
        encode_temp = polar_encoder(u, lambda_offset, llr_layer_vec);
    
        % bpsk modulation
        encode_temp = 1 - 2 * encode_temp;
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        llr = 2/sigma^2*receive_sample;

        decision_bits = CASCL_decoder(llr, L, K, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        
        count = sum(decision_bits' ~= source_crc_bit);
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end
        if (PerNum>=100 && iter>=10000)
            break;
        end
    end
    iterNum(i) = iter;
    per(i) = PerNum/iter;
    ber(i) = BerNum/(K*iter);
end