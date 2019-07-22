% clc
clear

% 基本参数设置
n = 8;  % 比特位数
N = 2^n;
R = 0.4531;    % 码率
Ng = 12;
poly = [1 1 1 1 1 0 0 0 1 0 0 1 1];

SNR = [4 4.5 5];
% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
Kp = 24;

K = 140;
Kinfo = 128;
Kpure = 104;
k_f = N-K;% frozen_bits length
G = encoding_matrix(n);


% get generate matrix
frozen_bits = zeros(1,k_f);

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
        
       
%         if mod(iter,100)==0
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum1: %2d\tNow PerNum2: %2d\tNow Error Bits: %2d', iter, SNR(i),PerNum1,PerNum2,BerNum1+BerNum2);
%         end
        
        
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
        
        decision_bits1 = polarSC_decoder(n,receive_sample1,sigma,frozen_index,frozen_bits,info_index);
        decision_bits2 = polarSC_decoder(n,receive_sample2,sigma,frozen_index,frozen_bits,info_index);
        
        
        err1 = sum(crccheck(decision_bits1,poly))==0;
        err2 = sum(crccheck(decision_bits2,poly))==0;
        % crc Check Result：If only one polar is uncorrect,then using BP
        % decoder with some concatenated bits extrasinc information.
        
        % situation 1: polar1 wrong, polr2 right;
        if ~err1 && err2
            % modify polar1 frozen_index frozen_bits info_index
            ReSC_oddWrong = ReSC_oddWrong + 1;
            %将K_p mutual bits 当作冻结位处理（直接附在原冻结位后面即可）
            frozen_index(k_f+1:k_f+Kp) = inter_index;
            frozen_bits(k_f+1:k_f+Kp) = mutual_bits;
            decision_bits1 = polarSC_decoder(n,receive_sample1,sigma,frozen_index,frozen_bits,info_index);
            if sum(crccheck(decision_bits1,poly)) == 0
               ReSC_oddCorrect =  ReSC_oddCorrect + 1;
               odd_correctFlag = true;  %设置flag,即如果此处CRC结果表示Redecoding正确了，但是后面统计PER出错,odd_missCheck++
            end
        end
        
        % situation 2: polar1 right, polr2 wrong;
        if err1 && ~err2
            ReSC_evenWrong = ReSC_evenWrong + 1;
            frozen_index(k_f+1:k_f+Kp) = inter_index;
            frozen_bits(k_f+1:k_f+Kp) = mutual_bits;
            decision_bits2 = polarSC_decoder(n,receive_sample2,sigma,frozen_index,frozen_bits,info_index);
            if sum(crccheck(decision_bits2,poly)) == 0
               ReSC_evenCorrect =  ReSC_evenCorrect + 1;
               even_correctFlag = true;
            end
        end
        
        if ~err1 && ~err2
            AllWrong = AllWrong + 1;
        end
        
        if err1 && err2
            allright_flag = true;
            AllRight = AllRight + 1;
        end
        
        % reset frozen_index and frozen_bits
        frozen_index = frozen_index(1:k_f);
        frozen_bits = frozen_bits(1:k_f);
        % situation 3 and 4: polar1 and polar2 are both right or wrong
        % we have no salution.
        
        % calculate BER and PER
        count1 = sum(decision_bits1(1:Kinfo) ~= source_bit1);
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
        count2 = sum(decision_bits2(1:Kinfo) ~= source_bit2);
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

path = './results/';
filename = [path, 'StandSC_N',num2str(N),'_R',num2str(R),'.mat'];
save(filename)