clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.5;    % 码率
Ng = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];

SNR = [0 1 2 3 3.5];
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
K = N*R;  % information bit length
Kp = N*R*0.25;  % Cascaded decoding length
k_f = N-K;% frozen_bits length
rng('shuffle');
for i = 1:length(SNR)
    
    %R = 0.25
%     filename = 'Pe_N256_snr3.2_R5.mat';
%     if SNR(i) == 3.5
%         filename = sprintf('Pe_N256_snr%1.1f_R5.mat',SNR(i));
%     else
%         filename = sprintf('Pe_N256_snr%d_R5.mat',SNR(i));
%     end
    
    % get information bits and concatenated bits
    filename = 'GA_N1024_R5_snr3.2.mat';
    % get information bits and concatenated bits
    load(filename);   % load the channel information
    [Ptmp, I] = sort(P,'descend');
    info_index = sort(I(1:K));  % 挑选质量好的信道传输信息位
    info_without_crc = info_index(1:K-Ng);  %得到K_{info}个信息位信道
    frozen_index = sort(I(K+1:end));   % 传输冻结位的信道
    inter_index = sort(I(K-Kp+1:K));
    
    % get generate matrix
    G = encoding_matrix(n);
    Gi = G(info_index,:);
    Gf = G(frozen_index,:);
    frozen_bits = randi([0 1],1,k_f);
    
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
    even_missCheck = 0;
    odd_missCheck = 0;
    allright_missCheckOdd = 0;
    allright_missCheckEven = 0;
    while true
        iter = iter +1;
        % 异常处理
        evenNum = (ReBP_evenWrong - ReBP_evenCorrect) + AllWrong;
        oddNum = (ReBP_oddWrong - ReBP_oddCorrect) + AllWrong;
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
        encode_temp1 = rem(source_crc_bit1*Gi + frozen_bits*Gf,2);
        encode_temp2 = rem(source_crc_bit2*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp1 = (-1).^(encode_temp1 + 1);
        encode_temp2 = (-1).^(encode_temp2 + 1);
        % add noise
        receive_sample1 = encode_temp1 + sigma * randn(size(encode_temp1));
        receive_sample2 = encode_temp2 + sigma * randn(size(encode_temp2));
        
        % BP decoder follow
        [lr_u1,lr_x1] = getBP_Parameter(receive_sample1,frozen_bits,frozen_index,n,init_max,sigma);
        decision_bits1 = polarBP_decoder(n,lr_u1,lr_x1,max_iter,info_index);
        [lr_u2,lr_x2] = getBP_Parameter(receive_sample2,frozen_bits,frozen_index,n,init_max,sigma);
        decision_bits2 = polarBP_decoder(n,lr_u2,lr_x2,max_iter,info_index);
        
        % CRC check follow
        receive_crc_bits1 = crccheck(decision_bits1,poly);
        receive_crc_bits2 = crccheck(decision_bits2,poly);
        % crc Check Result：If only one polar is uncorrect,then using BP
        % decoder with some concatenated bits extrasinc information.
        
        % reset the flag
        odd_correctFlag = false;
        even_correctFlag = false;
        allright_flag = false;
        
        % situation 1: polar1 wrong, polr2 right;
        if ~isempty(find(receive_crc_bits1,1)) && isempty(find(receive_crc_bits2,1))
            ReBP_oddWrong = ReBP_oddWrong + 1;
            for m = 1:length(temp_index)
                if decision_bits2(temp_index(m)) == 0
                    lr_u1(reverse_index(n,info_without_crc(temp_index(m)))) = init_max;
                else
                    lr_u1(reverse_index(n,info_without_crc(temp_index(m)))) = -init_max;
                end
            end
            decision_bits1 = polarBP_decoder(n,lr_u1,lr_x1,max_iter,info_index);
            if sum(crccheck(decision_bits1,poly)) == 0
               ReBP_oddCorrect =  ReBP_oddCorrect + 1;
               odd_correctFlag = true;
            end
        end
        
        % situation 2: polar1 right, polr2 wrong;
        if isempty(find(receive_crc_bits1,1)) && ~isempty(find(receive_crc_bits2,1))
            ReBP_evenWrong = ReBP_evenWrong + 1;
           for m = 1:length(temp_index)
                if decision_bits1(temp_index(m)) == 0
                    lr_u2(reverse_index(n,info_without_crc(temp_index(m)))) = init_max;
                else
                    lr_u2(reverse_index(n,info_without_crc(temp_index(m)))) = -init_max;
                end
            end
            decision_bits2 = polarBP_decoder(n,lr_u2,lr_x2,max_iter,info_index);
            if sum(crccheck(decision_bits2,poly)) == 0
               ReBP_evenCorrect =  ReBP_evenCorrect + 1;
               even_correctFlag = true;
            end
        end
        
        % situation 3 and 4: polar1 and polar2 are both right or wrong
        if ~isempty(find(receive_crc_bits1,1)) && ~isempty(find(receive_crc_bits2,1))
            AllWrong = AllWrong + 1;
        end
        
        if isempty(find(receive_crc_bits1,1)) && isempty(find(receive_crc_bits2,1))
            AllRight = AllRight + 1;
            allright_flag = true;
        end
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
    rs_oddwrong(i) = ReBP_oddWrong;
    rs_evenwrong(i) = ReBP_evenWrong;
    rs_oddcorr(i) = ReBP_oddCorrect;
    rs_evencorr(i) = ReBP_evenCorrect;
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
end