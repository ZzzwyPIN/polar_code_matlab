clc
clear

% 基本参数设置
n = 10;  % 比特位数 
N = 2^n;
m = 3;
Ng = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
K = 553;
Kp = 74;
SNR = [2 2.5 2.8 3 3.2 3.5];
% 参数计算
R = floor((K-Ng-Kp/3)/N * 10)/10;
k_f = N-K;% frozen_bits length
snr = 10.^(SNR/10);
esn0 = snr * R;
N = 2^n;

lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

%channel sort
filename = 'Pe_N1024_snr2.mat';
load(filename);
[~, I] = sort(P);
pure_info_index = I(1:K-Kp-Ng);  % 挑选质量好的信道传输信息位
MUUB = I(K-Kp+1:K);  % Bit channel of the most unreliable unfrozen bits
crc_index = I(K-Kp-Ng+1:K-Kp); % Bit channel of the CRC bits.
frozen_index = I(K+1:end);   % 传输冻结位的信道
info_index = [MUUB pure_info_index crc_index];

G = spc_encoding(n);

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
    % 以下参数用来记录每个SNR点，论文中提到的case1-case4发生次数
    OneWrong = 0;
    TwoWrong = 0;
    ThreeWrong = 0;    %odd block have correct new rounds of SCL decoding
    OneCorrect = 0;   %even block have correct new rounds of SCL decoding
    TwoCorrect = 0;
    ThreeCorrect = 0;
    TwoBlockWrong = 0;
    AllRight = 0;
    AllWrong = 0;
 
    while true 
        
        iter = iter + 1;
        % reset the frozen bits and mutual bits
        frozen_bits = ones(N,1);
        mutual_bits = zeros(N,1);
        frozen_bits(info_index) = 0;
        

        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum1: %2d\tNow Error Bits: %2d',iter,SNR(i),PerNum1,BerNum1+BerNum2+BerNum3);

        
        source_bits = rand(m, K-Ng)>0.5;
        source_bits(end,1:Kp) = mod(sum(source_bits(1:end-1,1:Kp)), 2);
        
        %Add CRC
        source_bits_crc = crcadd_m(source_bits, poly);
        
        %按对应位置放置信息位
        u = zeros(m,N);
        u(:,info_index) = source_bits_crc;
        
        %encoding
        codewords = mod(u * G, 2);
        
        %BPSK
        tran_samples = 1 - 2 * codewords;
        
        %add noise
        receive_sample = tran_samples + sigma * randn(size(tran_samples));
        
        %Computing LLR
        LLR = 2/sigma^2 * receive_sample';

        %decoding
        decision_bits = zeros(K,m);
        decision_bits(:,1) = SC_decoder(LLR(:,1), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        decision_bits(:,2) = SC_decoder(LLR(:,2), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        decision_bits(:,3) = SC_decoder(LLR(:,3), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
    
        err1 = sum(crccheck_m(decision_bits(:,1)',poly))==0;
        err2 = sum(crccheck_m(decision_bits(:,2)',poly))==0;
        err3 = sum(crccheck_m(decision_bits(:,3)',poly))==0;
        
        % crc Check Result：
        flag = err1+err2+err3;
        
        switch flag
            case 3
                AllRight = AllRight + 1;
            case 2
                % if only the block one wrong:
                if ~err1
                    OneWrong = OneWrong + 1;
                    frozen_bits(MUUB) = 2;
                    mutual_bits(MUUB) = mod(decision_bits(1:Kp,2)+decision_bits(1:Kp,3),2);
                    decision_bits(:,1) = SC_decoder(LLR(:,1), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                    if sum(crccheck_m(decision_bits(:,1)',poly))==0
                       OneCorrect = OneCorrect + 1;
                    end
                elseif ~err2
                    TwoWrong = TwoWrong + 1;
                    frozen_bits(MUUB) = 2;
                    mutual_bits(MUUB) = mod(decision_bits(1:Kp,1)+decision_bits(1:Kp,3),2);
                    decision_bits(:,2) = SC_decoder(LLR(:,2), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                    if sum(crccheck_m(decision_bits(:,2)',poly))==0
                       TwoCorrect = TwoCorrect + 1;
                    end
                elseif ~err3
                    ThreeWrong = ThreeWrong + 1;
                    frozen_bits(MUUB) = 2;
                    mutual_bits(MUUB) = mod(decision_bits(1:Kp,1)+decision_bits(1:Kp,2),2);
                    decision_bits(:,3) = SC_decoder(LLR(:,3), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                    if sum(crccheck_m(decision_bits(:,3)',poly))==0
                       ThreeCorrect = ThreeCorrect + 1;
                    end
                end
            case 1
                TwoBlockWrong = TwoBlockWrong + 1;
            case 0
                AllWrong = AllWrong + 1;
        end
        
        % calculate BER and PER
        count1 = sum(decision_bits(:,1)' ~= source_bits_crc(1,:));
        if count1 ~= 0
            PerNum1 = PerNum1 + 1;
            BerNum1 = BerNum1 + count1;
        end
        
        count2 = sum(decision_bits(:,2)' ~= source_bits_crc(2,:));
        if count2 ~= 0
            PerNum2 = PerNum2 + 1;
            BerNum2 = BerNum2 + count2;
        end
        
        count3 = sum(decision_bits(:,3)' ~= source_bits_crc(3,:));
        if count3 ~= 0
            PerNum3 = PerNum3 + 1;
            BerNum3 = BerNum3 + count3;
        end
       
        
        if (PerNum1>=100 && PerNum2>=100 && PerNum3>=100 && iter>=10000)
            break;
        end
        if iter>=10000000
           break; 
        end
        
   end    
    iterNum(i) = iter;
    per(i) = (PerNum1+PerNum2+PerNum3)/(3*iter);
    ber(i) = (BerNum1+BerNum2+PerNum3)/(3*K-Kp)/iter;
    
    onewrong(i) = OneWrong;
    twowrong(i) = TwoWrong;
    threewrong(i) = ThreeWrong;
    
    onecorr(i) = OneCorrect;
    twocorr(i) = TwoCorrect;
    threecorr(i) = ThreeCorrect;
    
    twoblockwrong(i) = TwoBlockWrong;
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
    
end
% record simulation results
% path = './results/';
% filename = [path, 'PCM_FastSC_N',num2str(N),'_R',num2str(R),'.mat'];
% save(filename)