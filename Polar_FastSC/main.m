clc
clear

% 基本参数设置
n = 10;  % 比特位数
R = 0.5;    % 码率
Ng = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
% L = 8;   %SCL List

SNR = [0 1 2 3 4];
% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
N = 2^n;

lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

K = N*R;  % information bit length
Kp = N*R*0.25;  % Cascaded decoding length
k_f = N-K;% frozen_bits length


%CRC
% [gen, det] = get_crc_objective(Ng);
% source_block = 2*k-k1;
% frozen_block = 2*k_f;
filename = 'GA_N1024_R5_snr3.2.mat'; 
% get information bits and concatenated bits
load(filename);   % load the channel information
[Ptmp, I] = sort(P,'descend');
info_index = sort(I(1:K));  % 挑选质量好的信道传输信息位

rng('shuffle');
for i = 1:length(SNR)
    
    sigma = (2*esn0(i))^(-0.5);
    % set PER and BER counter
    PerNum = 0;
    BerNum = 0;
    iter = 0;
    %counter the number of Re-SC decoding
    % 以下参数用来记录每个SNR点，论文中提到的case1-case4发生次数
 
    while true 
        
        
        iter = iter + 1;
        % reset the frozen bits and mutual bits
        frozen_bits = ones(N,1);
        frozen_bits(info_index) = 0;
        
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum: %2d\tNow Error Bits: %2d',iter,SNR(i),PerNum,BerNum);
        source_bit = rand(1,K)>0.5;


        u = zeros(N, 1);
        u(info_index) = source_bit;
        encode_temp = polar_encoder(u, lambda_offset, llr_layer_vec);

    
        % bpsk modulation
        encode_temp = 1 - 2 * encode_temp;

        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        llr = 2/sigma^2*receive_sample;
        
        
        decision_bits = SC_decoder(llr, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);


        count = sum(decision_bits' ~= source_bit);
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end

        if (PerNum>=10000 && iter>=100)
            break;
        end
        
        
    end    
    iterNum(i) = iter;
    per(i) = (PerNum)/(2*iter);
    ber(i) = (BerNum)/K/iter;

    
end