clc
clear

% 基本参数设置
n = 11;  % 比特位数
Ng = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
L = 32;   %SCL List
K = 1060; %the number of information bits of the underlying blocks
Kp = 40; %the number of mutual bits
SNR = 0.2:0.2:2;

%Compute the parameters
N = 2^n;
R = (K-Ng-Kp/2)/N;
snr = 10.^(SNR/10);
esn0 = snr * R;
k_f = N-K;% frozen_bits length

%pre-computation of some parameters in the decoding
lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

% selection of the bit channel.
load('Pe_N2048_snr2.mat');
[~, I] = sort(P);
pure_info_index = I(1:K-Kp-Ng);  % 挑选质量好的信道传输信息位
MUUB = I(K-Kp+1:K);  % Bit channel of the most unreliable unfrozen bits
crc_index = I(K-Kp-Ng+1:K-Kp); % Bit channel of the CRC bits.
frozen_index = I(K+1:end);   % 传输冻结位的信道

info_index = [MUUB pure_info_index crc_index];

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
    ReSCL_oddWrong = 0;
    ReSCL_evenWrong = 0;
    ReSCL_oddCorrect = 0;    %odd block have correct new rounds of SCL decoding
    ReSCL_evenCorrect = 0;   %even block have correct new rounds of SCL decoding
    AllRight = 0;
    AllWrong = 0;
    
    while true
        
        iter = iter + 1;
        % reset the frozen bits and mutual bits
        frozen_bits = ones(N,1);
        mutual_bits = zeros(N,1);
        frozen_bits(info_index) = 0;
        

        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum1: %2d\tNow PerNum2: %2d\tNow Error Bits: %2d', iter, SNR(i),PerNum1,PerNum2,BerNum1+BerNum2);

        
        %Generation of the source bits
        source_bit1 = rand(1,K-Ng)>0.5;
        
        % Insert the mutual bits
        source_bit2 = zeros(size(source_bit1));
        source_bit2(1:Kp) = source_bit1(1:Kp);
        source_bit2(Kp+1:end) = rand(1,K-Kp-Ng)>0.5;
        
        % CRC attachment.
        source_crc_bit1 = crcadd_m(source_bit1,poly);
        source_crc_bit2 = crcadd_m(source_bit2,poly);
        
        u1 = zeros(N, 1);
        u2 = zeros(N, 1);
        u1(info_index) = source_crc_bit1;
        u2(info_index) = source_crc_bit2;
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

        [decision_bits1, err1] = CASCL_decoder(llr1, L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        [decision_bits2, err2] = CASCL_decoder(llr2, L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        
        if err1
            err1 = sum(crccheck_m(decision_bits1',poly));
        end
      
        if err2
            err2 = sum(crccheck_m(decision_bits2',poly));
        end
       
        % situation 1: polar1 wrong, polr2 right;
        if err1 && ~err2
            ReSCL_oddWrong = ReSCL_oddWrong + 1;
            % modify polar1 frozen_index frozen_bits info_index
            frozen_bits(MUUB) = 2;
            mutual_bits(MUUB) = decision_bits2(1:Kp);
            [decision_bits1, ~] = CASCL_decoder(llr1, L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
            if sum(source_crc_bit1 ~= decision_bits1') == 0
                ReSCL_oddCorrect = ReSCL_oddCorrect + 1;
            end
        end
        
        % situation 2: polar1 right, polr2 wrong;
        %In this case, no need for new round of decoding therein.
        if ~err1 && err2
            ReSCL_evenWrong = ReSCL_evenWrong + 1;
            % modify polar1 frozen_index frozen_bits info_index
            frozen_bits(MUUB) = 2;
            mutual_bits(MUUB) = decision_bits1(1:Kp);
            [decision_bits2, ~] = CASCL_decoder(llr2, L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
            if sum(source_crc_bit2 ~= decision_bits2') == 0
                ReSCL_evenCorrect = ReSCL_evenCorrect + 1;
            end
        end
        
        %situation 3: All wrong, we have no solution
        if err1 && err2
            AllWrong = AllWrong + 1;
        end
        
        %situation 4: All right, just output
        if ~err1 && ~err2
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
        
        
        if (PerNum1>=100 && PerNum2>=100 && iter>=10000)
            break;
        end
        
        if iter >= 10000000
           break; 
        end
        
    end
    iterNum(i) = iter;
    per(i) = (PerNum1+PerNum2)/(2*iter);
    ber(i) = (BerNum1+BerNum2)/(2*K-Kp)/iter;
    rs_oddwrong(i) = ReSCL_oddWrong;
    rs_evenwrong(i) = ReSCL_evenWrong;
    rs_oddcorr(i) = ReSCL_oddCorrect;
    rs_evencorr(i) = ReSCL_evenCorrect;
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
end

% record the results
path = '../results/';
filename = [path, 'PCM_SCL',num2str(L),'_N',num2str(N),'_R',num2str(R),'.mat'];
save(filename)