clc
clear
addpath('../GA/')
% 基本参数设置
n = 8;  % 比特位数
N = 2^n;
R = 0.4531;    % 码率
Ng = 12;
poly = [1 1 1 1 1 0 0 0 1 0 0 1 1];

SNR = [2 3 4 4.5 5];
init_lr_max = 3;    % limit the max LR of the channel to be with [-3 3]
max_iter = 40;

% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
init_max = init_lr_max * n;
if init_max > 30
    init_max = 30;
end
Kp = 24;

K = 140;
Kinfo = 128;
Kpure = 104;
k_f = N-K;% frozen_bits length
G = encoding_matrix(n);

rng('shuffle');
for i = 1:length(SNR)
    sigma = (2*esn0(i))^(-0.5);

    P = GA(sigma,N);
    [~, I] = sort(P,'descend');
    pure_index = I(1:Kpure);
    inter_index = I(117:140);
    crc_index = I(105:116);
    info_index = [pure_index inter_index crc_index];
    frozen_index = I(K+1:end);   % 传输冻结位的信道
    
    % get generate matrix
    
    
    
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
        
        if mod(iter,100) == 0
            fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum1: %2d\tNow PerNum2: %2d\tNow Error Bits: %2d', iter, SNR(i),PerNum1,PerNum2,BerNum1+BerNum2);
        end
        source_bit1 = randi([0 1],1,Kinfo);
        mutual_bits = source_bit1(Kpure+1:Kinfo);
        source_bit2 = randi([0 1],1,Kpure);
        source_bit2(Kpure+1:Kinfo) = mutual_bits;
        source_crc_bit1 = crcadd(source_bit1,poly);
        source_crc_bit2 = crcadd(source_bit2,poly);
        
        u1 = zeros(1,N);
        u2 = zeros(1,N);
        u1(info_index) = source_crc_bit1;
        u2(info_index) = source_crc_bit2;
        
        encode_bits1 = rem(u1*G,2);
        encode_bits2 = rem(u2*G,2);

        % bpsk modulation
        encode_temp1 = 1 - 2 * encode_bits1;
        encode_temp2 = 1 - 2 * encode_bits2;
        % add noise
        receive_sample1 = encode_temp1 + sigma * randn(size(encode_temp1));
        receive_sample2 = encode_temp2 + sigma * randn(size(encode_temp2));
        
        % BP decoder follow
        lr_x1 = 2*receive_sample1/(sigma^2);
        % decoding follow
        lr_u1 = zeros(1,N); % save send sample LR in each iteration
        lr_u1(reverse_index(n,frozen_index)) = init_max;
        decision_bits1 = polarBP_decoder(n,lr_u1,lr_x1,max_iter,info_index);
        
        % block even decoding
        lr_x2 = 2*receive_sample2/(sigma^2);
        % decoding follow
        lr_u2 = zeros(1,N); % save send sample LR in each iteration
        lr_u2(reverse_index(n,frozen_index)) = init_max;
        decision_bits2 = polarBP_decoder(n,lr_u2,lr_x2,max_iter,info_index);
        
        % CRC check follow
        err1 = sum(crccheck(decision_bits1,poly))==0;
        err2 = sum(crccheck(decision_bits2,poly))==0;
        % crc Check Result：If only one polar is uncorrect,then using BP
        % decoder with some concatenated bits extrasinc information.
        
        % reset the flag
        odd_correctFlag = false;
        even_correctFlag = false;
        allright_flag = false;
        
        % situation 1: polar1 wrong, polr2 right;
        if ~err1 && err2
            ReBP_oddWrong = ReBP_oddWrong + 1;
            for m = 1:Kp   
                lr_u1(reverse_index(n,inter_index(m))) = (1-2*mutual_bits(m))*init_max;  
            end
            decision_bits1 = polarBP_decoder(n,lr_u1,lr_x1,max_iter,info_index);
            if sum(decision_bits1~=source_crc_bit1) == 0
               ReBP_oddCorrect =  ReBP_oddCorrect + 1;
               odd_correctFlag = true;
            end
        end
        
        % situation 2: polar1 right, polr2 wrong;
        if err1 && ~err2
            ReBP_evenWrong = ReBP_evenWrong + 1;
            for m = 1:Kp   
                lr_u2(reverse_index(n,inter_index(m))) = (1-2*mutual_bits(m))*init_max;  
            end
            decision_bits2 = polarBP_decoder(n,lr_u2,lr_x2,max_iter,info_index);
            if sum(decision_bits2~=source_crc_bit2) == 0
               ReBP_evenCorrect =  ReBP_evenCorrect + 1;
               even_correctFlag = true;
            end
        end
        
        % situation 3 and 4: polar1 and polar2 are both right or wrong
        if ~err1 && ~err2
            AllWrong = AllWrong + 1;
        end
        
        if err1 && err2
            AllRight = AllRight + 1;
            allright_flag = true;
        end
        % we have no salution.
        
        % calculate BER and PER
        count1 = sum(decision_bits1(1:Kinfo)~=source_bit1);
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
        count2 = sum(decision_bits2(1:Kinfo)~=source_bit2);
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
    ber(i) = (BerNum1+BerNum2)/(2*Kinfo-Kp)/iter;
    rs_oddwrong(i) = ReBP_oddWrong;
    rs_evenwrong(i) = ReBP_evenWrong;
    rs_oddcorr(i) = ReBP_oddCorrect;
    rs_evencorr(i) = ReBP_evenCorrect;
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
end

path = './results/';
filename = [path, 'StandBP_N',num2str(N),'_R',num2str(R),'.mat'];
save(filename)