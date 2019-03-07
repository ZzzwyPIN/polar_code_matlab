clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.5;    % 码率
Ng = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];

SNR = [0 1 2 3 3.2];
% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
N = 2^n;
K = N*R;  % information bit length
Kp = N*R*0.25;  % Cascaded decoding length
k_f = N-K;% frozen_bits length
% source_block = 2*k-k1;
% frozen_block = 2*k_f;
filename = 'Pe_N256_snr3.2_R5.mat'; 
% get information bits and concatenated bits
load(filename);   % load the channel information
[Ptmp, I] = sort(P);
info_index = sort(I(1:K));  % 挑选质量好的信道传输信息位
info_without_crc = info_index(1:K-Ng);
frozen_index = sort(I(K+1:end));   % 传输冻结位的信道
[~,temp] = sort(P(info_without_crc));
inter_index = sort(info_without_crc(temp(end:-1:end-Kp+1)));
clear temp;

% get generate matrix
G = encoding_matrix(n);
Gi = G(info_index,:);
Gf = G(frozen_index,:);
frozen_bits = randi([0 1],1,k_f);

rng('shuffle');
for i = 1:length(SNR)
    sigma = (2*esn0(i))^(-0.5);
    % set PER and BER counter
    PerNum1 = 0;
    BerNum1 = 0;
    PerNum2 = 0;
    BerNum2 = 0;
    PerNum3 = 0;
    BerNum3 = 0;
    iter = 0;
    %counter the number of Re-SC decoding
    Block1_ReCorrect = 0;  %New rounds of SC decoding correctly
    Block2_ReCorrect = 0;
    Block3_ReCorrect = 0;
    Block1_Wrong = 0;  %The number of only one block decoding incorrect
    Block2_Wrong = 0;
    Block3_Wrong = 0;
    TwoBlockWrong = 0;
    
    AllRight = 0;
    AllWrong = 0;
    while true
        
        iter = iter + 1;
        fprintf('\nNow iter: %2d\tSNR: %d\tPerNum1: %2d\tPerNum2: %2d\tPerNum3: %2d', iter,SNR(i),PerNum1,PerNum2,PerNum3);
        source_bit1 = randi([0 1],1,K-Ng);
        source_bit2 = randi([0 1],1,K-Ng);
        source_bit3 = randi([0 1],1,K-Kp-Ng);
        
        [~,temp_index] = ismember(inter_index,info_without_crc);
        mutualBits = rem(source_bit1(temp_index)+source_bit2(temp_index),2);
        
        source_bit3 = insert_bit(mutualBits,source_bit3,temp_index);
        source_crc_bit1 = crcadd(source_bit1,poly);
        source_crc_bit2 = crcadd(source_bit2,poly);
        source_crc_bit3 = crcadd(source_bit3,poly);
        encode_temp1 = rem(source_crc_bit1*Gi + frozen_bits*Gf,2);
        encode_temp2 = rem(source_crc_bit2*Gi + frozen_bits*Gf,2);
        encode_temp3 = rem(source_crc_bit3*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp1 = (-1).^(encode_temp1 + 1);
        encode_temp2 = (-1).^(encode_temp2 + 1);
        encode_temp3 = (-1).^(encode_temp3 + 1);
        % add noise
        receive_sample1 = encode_temp1 + sigma * randn(size(encode_temp1));
        receive_sample2 = encode_temp2 + sigma * randn(size(encode_temp2));
        receive_sample3 = encode_temp3 + sigma * randn(size(encode_temp3));
        
        decision_bits1 = polarSC_decoder(n,receive_sample1,sigma,frozen_index,frozen_bits,info_index);
        decision_bits2 = polarSC_decoder(n,receive_sample2,sigma,frozen_index,frozen_bits,info_index);
        decision_bits3 = polarSC_decoder(n,receive_sample3,sigma,frozen_index,frozen_bits,info_index);
        
        
        receive_crc_bits1 = crccheck(decision_bits1,poly);
        receive_crc_bits2 = crccheck(decision_bits2,poly);
        receive_crc_bits3 = crccheck(decision_bits3,poly);
        % crc Check Result：
        % decoder with some concatenated bits extrasinc information.
%         newRoundFlag = false;


        block_flag = [~isempty(find(receive_crc_bits1,1)) ~isempty(find(receive_crc_bits2,1)) ~isempty(find(receive_crc_bits3,1))];
        % situation: Only one block decoding incorrectly.
        switch sum(block_flag)
            case 0
                AllRight = AllRight + 1;
            case 1
                if find(block_flag) == 1
                    Block1_Wrong = Block1_Wrong + 1;
                    mutualVector = rem(decision_bits2+decision_bits3,2);
                    [frozen_index,frozen_bits] = modifyFrozenIndexAndBits(frozen_index,frozen_bits,info_without_crc,temp_index,mutualVector);
                    decision_bits1 = polarSC_decoder(n,receive_sample1,sigma,frozen_index,frozen_bits,info_index);
                    if sum(crccheck(decision_bits1,poly))==0
                        Block1_ReCorrect = Block1_ReCorrect+1;
                    end
                elseif find(block_flag)==2
                    Block2_Wrong = Block2_Wrong + 1;
                    mutualVector = rem(decision_bits1+decision_bits3,2);
                    [frozen_index,frozen_bits] = modifyFrozenIndexAndBits(frozen_index,frozen_bits,info_without_crc,temp_index,mutualVector);
                    decision_bits2 = polarSC_decoder(n,receive_sample2,sigma,frozen_index,frozen_bits,info_index);
                    if sum(crccheck(decision_bits2,poly))==0
                        Block2_ReCorrect = Block2_ReCorrect+1;
                    end
                elseif find(block_flag)==3
                    Block3_Wrong = Block3_Wrong + 1;
                    mutualVector = rem(decision_bits1+decision_bits2,2);
                    [frozen_index,frozen_bits] = modifyFrozenIndexAndBits(frozen_index,frozen_bits,info_without_crc,temp_index,mutualVector);
                    decision_bits3 = polarSC_decoder(n,receive_sample3,sigma,frozen_index,frozen_bits,info_index);
                    if sum(crccheck(decision_bits3,poly))==0
                        Block3_ReCorrect = Block3_ReCorrect+1;
                    end
                end
            case 2
                %we dont do anything
                TwoBlockWrong = TwoBlockWrong + 1;
            case 3
                AllWrong = AllWrong + 1; 
        end
        
        
        frozen_index = frozen_index(1:k_f);
        frozen_bits = frozen_bits(1:k_f);
        % situation 3 and 4: polar1 and polar2 are both right or wrong
        % we have no salution.
        
        % calculate BER and PER
        counter = [sum(decision_bits1 ~= source_crc_bit1) sum(decision_bits2 ~= source_crc_bit2) sum(decision_bits3 ~= source_crc_bit3)];
    
        
        
        if counter(1) ~= 0
            PerNum1 = PerNum1 + 1;
            BerNum1 = BerNum1 + counter(1);
            %CRC missCheck counter
        end
        
        if counter(2) ~= 0
            PerNum2 = PerNum2 + 1;
            BerNum2 = BerNum2 + counter(2);
        end
        
        if counter(3) ~= 0
            PerNum3 = PerNum3 + 1;
            BerNum3 = BerNum3 + counter(3);
        end
        
        if (PerNum1>=100 && PerNum2>=100 && PerNum3>=100 && iter>=10000)
            break;
        end
        
        
        finalWrong1 = Block1_Wrong-Block1_ReCorrect;
        finalWrong2 = Block2_Wrong-Block2_ReCorrect;
        finalWrong3 = Block3_Wrong-Block3_ReCorrect;
        
        if (3*AllWrong+2*TwoBlockWrong+finalWrong1+finalWrong2+finalWrong3)~=(PerNum1+PerNum2+PerNum3)
            warning('counter error:CRC check may not correct!');
            break;
        end
        
    end
    iterNum(i) = iter;
    per(i) = (PerNum1+PerNum2+PerNum3)/(3*iter);
    ber(i) = (BerNum1+BerNum2+BerNum3)/(3*K-Kp)/iter;
    
    block1_wrong(i) = Block1_Wrong;
    block1_reco(i) = Block1_ReCorrect;
    block2_wrong(i) = Block2_Wrong;
    block2_reco(i) = Block2_ReCorrect;
    block3_wrong(i) = Block3_Wrong;
    block3_reco(i) = Block3_ReCorrect;
    
    
    two_block_wrong(i) = TwoBlockWrong+1;
    
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
end