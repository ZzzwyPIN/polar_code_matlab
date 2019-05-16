clc
clear

% 基本参数设置
n = 10;  % 比特位数
R = 0.5;    % 码率

Ng = 16; %CRC bit length
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];

SNR = [1 2 3 4];
max_iter = 40;

% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;

N = 2^n;
K = floor(N*R);  % information bit length
k_f = N - K;

Kp = floor(K*0.25);  % Cascaded decoding length

[M_up, M_down] = index_Matrix(N);
lambda_offset = 2.^(0 : n);
llr_layer_vec = get_llr_layer(N);

% get information bits and concatenated bits
load('GA_N1024_R5_snr3.2.mat');   % load the channel information
[Ptmp, I] = sort(P,'descend');
info_index = sort(I(1:K));  % 挑选质量好的信道传输信息位
info_without_crc = info_index(1:K-Ng);  %得到K_{info}个信息位信道
frozen_index = sort(I(K+1:end));   % 传输冻结位的信道
inter_index = sort(I(K-Kp+1:K));

rng('shuffle')
for i = 1:length(SNR)
    sigma = (2*esn0(i))^(-0.5);
    % set PER and BER counter
    PerNum1 = 0;
    BerNum1 = 0;
    PerNum2 = 0;
    BerNum2 = 0;
    
    ReBP_oddWrong = 0;
    ReBP_evenWrong = 0;
    ReBP_oddCorrect = 0;
    ReBP_evenCorrect = 0;
    AllWrong = 0;
    AllRight = 0;
    
    iter = 0;
    while true
        
        
        iter = iter + 1;
        
        frozen_bits = ones(N , 1);
        frozen_bits(info_index) = 0;
        mutual_bits = zeros(N,1);
        
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow perNum1: %2d\tNow perNum2: %2d\tNow BerNum: %2d', iter, SNR(i),PerNum1,PerNum2,BerNum1+BerNum2);
        u1 = zeros(N, 1);
        u2 = zeros(N, 1);
        source_bit1 = randi([0 1],1,K-Ng);
        source_bit2 = randi([0 1],1,K-Kp-Ng);
        [~,temp_index] = ismember(inter_index,info_without_crc);
        source_bit2 = insert_bit(source_bit1,source_bit2,temp_index,temp_index);
        source_crc_bit1 = crcadd(source_bit1,poly);
        source_crc_bit2 = crcadd(source_bit2,poly);
        
        
        u1(info_index) = source_crc_bit1;
        u2(info_index) = source_crc_bit2;
        encode_temp1 = polar_encoder(u1, lambda_offset, llr_layer_vec);
        encode_temp2 = polar_encoder(u1, lambda_offset, llr_layer_vec);
    
        % bpsk modulation
        encode_temp1 = 1 - 2*encode_temp1;
        encode_temp2 = 1 - 2*encode_temp2;
        % add noise
        receive_sample1 = encode_temp1 + sigma * randn(size(encode_temp1));
        receive_sample2 = encode_temp2 + sigma * randn(size(encode_temp2));
        
        llr1 = 2*receive_sample1/(sigma^2);
        llr2 = 2*receive_sample2/(sigma^2);
        
        [decision_bits1, ~, ~, ~] = BP_Decoder_LLR(info_index, frozen_bits, llr1, max_iter, M_up, M_down);
        [decision_bits2, ~, ~, ~] = BP_Decoder_LLR(info_index, frozen_bits, llr2, max_iter, M_up, M_down);
        
        err1 = sum(crccheck(decision_bits1,poly))==0;
        err2 = sum(crccheck(decision_bits2,poly))==0;
        
        if ~err1 && err2
            ReBP_oddWrong = ReBP_oddWrong + 1;
            frozen_bits(inter_index) = 2;
            mutual_bits(inter_index) = decision_bits2(temp_index);
            [decision_bits1, ~, ~, ~] = BP_Decoder_LLR(info_index, frozen_bits, llr1, max_iter, M_up, M_down,mutual_bits);
            if sum(source_crc_bit1 ~= decision_bits1') == 0
                ReBP_oddCorrect = ReBP_oddCorrect + 1;
            end
        end
        
        if err1 && ~err2
            ReBP_evenWrong = ReBP_evenWrong + 1;
            frozen_bits(inter_index) = 2;
            mutual_bits(inter_index) = decision_bits1(temp_index);
            [decision_bits2, ~, ~, ~] = BP_Decoder_LLR(info_index, frozen_bits, llr2, max_iter, M_up, M_down,mutual_bits);
            if sum(source_crc_bit2 ~= decision_bits2') == 0
                ReBP_evenCorrect = ReBP_evenCorrect + 1;
            end
        end
        
        if ~err1 && ~err2
            AllWrong = AllWrong + 1;
        end
        
        if err1 && err2
            allright_flag = true;
            AllRight = AllRight + 1;
        end
        
        % calculate BER and PER
        count1 = sum(decision_bits1' ~= source_crc_bit1);
        if count1 ~= 0
            PerNum1 = PerNum1 + 1;
            BerNum1 = BerNum1 + count1;
        end
        count2 = sum(decision_bits2' ~= source_crc_bit2);
        if count2 ~= 0
            PerNum2 = PerNum2 + 1;
            BerNum2 = BerNum2 + count2;
        end
        
        if (PerNum1>=100 && PerNum2>=100 && iter>=1000)
            break;
        end
    end
    iterNum(i) = iter;
    per(i) = (PerNum1+PerNum2)/(2*iter);
    ber(i) = (BerNum1+BerNum2)/(2*K-Kp)/iter;
    rs_oddwrong(i) = ReBP_oddWrong;
    rs_evenwrong(i) = ReBP_evenWrong;
    rs_oddcorr(i) = ReBP_oddCorrect;
    rs_evencorr(i) = ReBP_evenCorrect;
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
end