% clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.5;    % 码率
Ng = 16;
L = 8;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];

SNR = 3;
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
info_without_crc = info_index(1:K-Ng);  %得到K_{info}个信息位信道
frozen_index = sort(I(K+1:end));   % 传输冻结位的信道

% [~,temp] = sort(P(info_without_crc));   %对K_{info}个信道再进行一次排序
inter_index = sort(I(K-Kp+1:K));
% inter_index = sort(info_without_crc(temp(end:-1:end-Kp+1)));    %取K_{info}中最差的K_p个信道
% clear temp;

lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

% get generate matrix
frozen_bits = zeros(1,k_f);
frozen_indices = ones(N,1);
frozen_indices(info_index) = 0;
info_bits_logical = logical(mod(frozen_indices + 1, 2));

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
    ReSC_oddCorrect = 0;    %odd block have correct new rounds of SC decoding
    ReSC_evenCorrect = 0;   %even block have correct new rounds of SC decoding
    AllRight = 0;
    AllWrong = 0;
    % 以下参数用下记录CRC校验出错的次数，即：CRC校验是正确的但是后面误组率统计却出错。
    even_missCheck = 0;
    odd_missCheck = 0;
    allright_missCheckOdd = 0;
    allright_missCheckEven = 0;
    
    while true 
        odd_correctFlag = false;
        even_correctFlag = false;
        allright_flag = false;
        
        iter = iter + 1;
        
        % 当出现CRC校验错误或者统计错误时直接跳出
        evenNum = (ReSC_evenWrong - ReSC_evenCorrect) + AllWrong;
        oddNum = (ReSC_oddWrong - ReSC_oddCorrect) + AllWrong;
        if (oddNum ~= PerNum1 || evenNum ~= PerNum2)
           break; 
        end
        
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum1: %2d\tNow oddNum: %2d\tNow PerNum2: %2d\tNow evenNum: %2d\tNow Error Bits: %2d', iter, SNR(i),PerNum1,oddNum,PerNum2,evenNum,BerNum1+BerNum2);
        source_bit1 = randi([0 1],1,K-Ng);
        source_bit2 = randi([0 1],1,K-Kp-Ng);
        [~,temp_index] = ismember(inter_index,info_without_crc);
        source_bit2 = insert_bit(source_bit1,source_bit2,temp_index,temp_index);
        source_crc_bit1 = crcadd(source_bit1,poly);
        source_crc_bit2 = crcadd(source_bit2,poly);
        u = zeros(N, 1);
        u(info_bits_logical) = source_crc_bit1;
        encode_temp1 = polar_encoder(u, lambda_offset, llr_layer_vec);
    
        % bpsk modulation
        encode_temp1 = 1 - 2 * encode_temp1;
        receive_sample1 = encode_temp1 + sigma * randn(size(encode_temp1));
        llr = 2/sigma^2 * receive_sample1';
        decision_bits = CASCL_decoder(llr, L, K, frozen_indices, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        
        
        
        
        
        
        
        
        
        
        
        encode_temp2 = (-1).^(encode_temp2 + 1);
        % add noise
        receive_sample1 = encode_temp1 + sigma * randn(size(encode_temp1));
        receive_sample2 = encode_temp2 + sigma * randn(size(encode_temp2));
        
        decision_bits1 = polarSC_decoder(n,receive_sample1,sigma,frozen_index,frozen_bits,info_index);
        decision_bits2 = polarSC_decoder(n,receive_sample2,sigma,frozen_index,frozen_bits,info_index);
        
        
        receive_crc_bits1 = crccheck(decision_bits1,poly);
        receive_crc_bits2 = crccheck(decision_bits2,poly);
        % crc Check Result：If only one polar is uncorrect,then using BP
        % decoder with some concatenated bits extrasinc information.
        
        % situation 1: polar1 wrong, polr2 right;
        if ~isempty(find(receive_crc_bits1,1)) && isempty(find(receive_crc_bits2,1))
            % modify polar1 frozen_index frozen_bits info_index
            ReSC_oddWrong = ReSC_oddWrong + 1;
            %将K_p mutual bits 当作冻结位处理（直接附在原冻结位后面即可）
            [frozen_index,frozen_bits] = modifyFrozenIndexAndBits(frozen_index,frozen_bits,info_without_crc,temp_index,decision_bits2);
            decision_bits1 = polarSC_decoder(n,receive_sample1,sigma,frozen_index,frozen_bits,info_index);
            if sum(crccheck(decision_bits1,poly)) == 0
               ReSC_oddCorrect =  ReSC_oddCorrect + 1;
               odd_correctFlag = true;  %设置flag,即如果此处CRC结果表示Redecoding正确了，但是后面统计PER出错,odd_missCheck++
            end
        end
        
        % situation 2: polar1 right, polr2 wrong;
        if isempty(find(receive_crc_bits1,1)) && ~isempty(find(receive_crc_bits2,1))
            ReSC_evenWrong = ReSC_evenWrong + 1;
            [frozen_index,frozen_bits] = modifyFrozenIndexAndBits(frozen_index,frozen_bits,info_without_crc,temp_index,decision_bits1);
            decision_bits2 = polarSC_decoder(n,receive_sample2,sigma,frozen_index,frozen_bits,info_index);
            if sum(crccheck(decision_bits2,poly)) == 0
               ReSC_evenCorrect =  ReSC_evenCorrect + 1;
               even_correctFlag = true;
            end
        end
        
        if ~isempty(find(receive_crc_bits1,1)) && ~isempty(find(receive_crc_bits2,1))
            AllWrong = AllWrong + 1;
        end
        
        if isempty(find(receive_crc_bits1,1)) && isempty(find(receive_crc_bits2,1))
            allright_flag = true;
            AllRight = AllRight + 1;
        end
        
        % reset frozen_index and frozen_bits
        frozen_index = frozen_index(1:k_f);
        frozen_bits = frozen_bits(1:k_f);
        % situation 3 and 4: polar1 and polar2 are both right or wrong
        % we have no salution.
        
        % calculate BER and PER
        count1 = sum(decision_bits1 ~= source_crc_bit1);
        if count1 ~= 0
            PerNum1 = PerNum1 + 1;
            BerNum1 = BerNum1 + count1;
            if odd_correctFlag
                odd_missCheck = odd_missCheck + 1;
            end
            if allright_flag
               allright_missCheckOdd = allright_missCheckOdd + 1; 
            end
        end
        count2 = sum(decision_bits2 ~= source_crc_bit2);
        if count2 ~= 0
            PerNum2 = PerNum2 + 1;
            BerNum2 = BerNum2 + count2;
            if even_correctFlag
                even_missCheck = even_missCheck + 1;
            end
            if allright_flag
               allright_missCheckEven = allright_missCheckEven + 1; 
            end
        end
        
        
        if (PerNum1>=100 && PerNum2>=100 && iter>=10000)
            break;
        end
        
        
        
    end
    iterNum(i) = iter;
    per(i) = (PerNum1+PerNum2)/(2*iter);
    ber(i) = (BerNum1+BerNum2)/(2*K-Kp)/iter;
    rs_oddwrong(i) = ReSC_oddWrong;
    rs_evenwrong(i) = ReSC_evenWrong;
    rs_oddcorr(i) = ReSC_oddCorrect;
    rs_evencorr(i) = ReSC_evenCorrect;
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
    odd_missCheckNum(i) = odd_missCheck;
    even_missCheckNum(i) = even_missCheck;
    allright_missCheckOddNum(i) = allright_missCheckOdd;
    allright_missCheckEvenNum(i) = allright_missCheckEven;
end