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
info_without_crc = info_index(1:K-Ng);  %得到K_{info}个信息位信道
frozen_index = sort(I(K+1:end));   % 传输冻结位的信道
inter_index = sort(I(K-Kp+1:K));

%get CRC check
% [~, ~, g] = get_crc_objective(Ng);
% [G_crc, H_crc] = crc_generator_matrix(g, K - Ng);
% crc_parity_check = G_crc(:, K - Ng + 1 : end)';
rng('shuffle');
for i = 1:length(SNR)
    
    sigma = (2*esn0(i))^(-0.5);
    % set PER and BER counter
    PerNum1 = 0;
    BerNum1 = 0;
    PerNum2 = 0;
    BerNum2 = 0;
    iter = 0;
    %counter the number of Re-SC decoding
    % 以下参数用来记录每个SNR点，论文中提到的case1-case4发生次数
    ReSC_oddWrong = 0;
    ReSC_evenWrong = 0;
    ReSC_oddCorrect = 0;    %odd block have correct new rounds of SCL decoding
    ReSC_evenCorrect = 0;   %even block have correct new rounds of SCL decoding
    AllRight = 0;
    AllWrong = 0;
 
    while true 
        
        
        iter = iter + 1;
        % reset the frozen bits and mutual bits
        frozen_bits = ones(N,1);
        mutual_bits = zeros(N,1);
        frozen_bits(info_index) = 0;
        info_bits_logical = logical(mod(frozen_bits + 1, 2));
        
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum1: %2d\tNow PerNum2: %2d\tNow Error Bits: %2d',iter,SNR(i),PerNum1,PerNum2,BerNum1+BerNum2);
        source_bit1 = rand(1,K-Ng)>0.5;
        source_bit2 = rand(1,K-Kp-Ng)>0.5;
        [~,temp_index] = ismember(inter_index,info_without_crc);
        source_bit2 = insert_bit(source_bit1,source_bit2,temp_index,temp_index);
        %info_with_crc = [info; mod(crc_parity_check * info, 2)];
        source_crc_bit1 = crcadd(source_bit1,poly);
        source_crc_bit2 = crcadd(source_bit2,poly);
%         source_crc_bit1 = [source_bit1; mod(crc_parity_check * source_bit1, 2)];
%         source_crc_bit2 = [source_bit2; mod(crc_parity_check * source_bit2, 2)];
        u1 = zeros(N, 1);
        u2 = zeros(N, 1);
        u1(info_bits_logical) = source_crc_bit1;
        u2(info_bits_logical) = source_crc_bit2;
        encode_temp1 = polar_encoder(u1, lambda_offset, llr_layer_vec);
        encode_temp2 = polar_encoder(u2, lambda_offset, llr_layer_vec);
    
        % bpsk modulation
        encode_temp1 = 1 - 2 * encode_temp1;
        encode_temp2 = 1 - 2 * encode_temp2;
        % add noise
        receive_sample1 = encode_temp1 + sigma * randn(size(encode_temp1));
        receive_sample2 = encode_temp2 + sigma * randn(size(encode_temp2));
        
        llr1 = 2/sigma^2*receive_sample1;
        llr2 = 2/sigma^2*receive_sample2;
        
        
        decision_bits1 = SC_decoder(llr1, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        decision_bits2 = SC_decoder(llr2, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);

%         [decision_bits1, err1] = CASCL_decoder(llr1, L, K, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
%         [decision_bits2, err2] = CASCL_decoder(llr2, L, K, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        
        err1 = sum(crccheck(decision_bits1,poly))==0;
        err2 = sum(crccheck(decision_bits2,poly))==0;
        % crc Check Result：If only one polar is uncorrect,then using BP
        % decoder with some concatenated bits extrasinc information.
        
        % situation 1: polar1 wrong, polr2 right;
        if ~err1 && err2
            % modify polar1 frozen_index frozen_bits info_index
            ReSC_oddWrong = ReSC_oddWrong + 1;
            % modify polar1 frozen_index frozen_bits info_index
            frozen_bits(inter_index) = 2;
            mutual_bits(inter_index) = decision_bits2(temp_index);
            decision_bits1 = SC_decoder(llr1, L, K, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
            if sum(source_crc_bit1 ~= decision_bits1') == 0
                ReSC_oddCorrect = ReSC_oddCorrect + 1;
            end
        end
        
        % situation 2: polar1 right, polr2 wrong;
        if err1 && ~err2
            ReSC_evenWrong = ReSC_evenWrong + 1;
            % modify polar1 frozen_index frozen_bits info_index
            frozen_bits(inter_index) = 2;
            mutual_bits(inter_index) = decision_bits1(temp_index);
            decision_bits2 = SC_decoder(llr2, L, K, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
            if sum(source_crc_bit2 ~= decision_bits2') == 0
                ReSC_evenCorrect = ReSC_evenCorrect + 1;
            end
        end
        
        if ~err1 && ~err2
            AllWrong = AllWrong + 1;
        end
        
        if err1 && err2
            allright_flag = true;
            AllRight = AllRight + 1;
        end
        
        
        % situation 3 and 4: polar1 and polar2 are both right or wrong
        % we have no salution.
        
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
        
        
        if (PerNum1>=10000 && PerNum2>=10000 && iter>=100)
            break;
        end
        
        
   end    
    enditerNum(i) = iter;
    per(i) = (PerNum1+PerNum2)/(2*iter);
    ber(i) = (BerNum1+BerNum2)/(2*K-Kp)/iter;
    rs_oddwrong(i) = ReSC_oddWrong;
    rs_evenwrong(i) = ReSC_evenWrong;
    rs_oddcorr(i) = ReSC_oddCorrect;
    rs_evencorr(i) = ReSC_evenCorrect;
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
    
end