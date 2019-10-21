%Improved PCM using SPC ecoding and decoding. And the punctured is used to
%keep the code rate the same as that of the stand-alone SC.

clc
clear

% 基本参数设置
n = 10;  % 比特位数
Ng = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
L = 8;   %SCL List
K = 548; %the number of information bits of the underlying blocks
Kp = 40; %the number of mutual bits
block_num = 2;
init_max = 3;
SNR = 0.5:0.5:3;

%Compute the parameters
N = 2^n;
R = (K-Ng-Kp/2)/N;
snr = 10.^(SNR/10);
esn0 = snr * R;
k_f = N-K;% frozen_bits length
init_lr_max = init_max * n;

%pre-computation of some parameters in the decoding
lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

% selection of the bit channel.
load('Pe_N1024_snr2.mat');
[~, I] = sort(P);
pure_info_index = I(1:K-Kp-Ng);  % 挑选质量好的信道传输信息位
MUUB = I(K-Kp+1:K);  % Bit channel of the most unreliable unfrozen bits
crc_index = I(K-Kp-Ng+1:K-Kp); % Bit channel of the CRC bits.
frozen_index = I(K+1:end);   % 传输冻结位的信道
info_index = [MUUB pure_info_index crc_index];

%Generation matrix
G = spc_encoding(n);
Gaa = G(info_index, info_index);
Gab = G(info_index, frozen_index);

%Set the frozen bits
frozen_bits = ones(N,1);
frozen_bits(info_index) = 0;

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
        
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum1: %2d\tNow PerNum2: %2d\tNow Error Bits: %2d', iter, SNR(i),PerNum1,PerNum2,BerNum1+BerNum2);
        
        %Generation of the source bits
        codeword = zeros(N,block_num);
        decision_bits = zeros(block_num,K);
        Xa = rand(block_num,K-Ng)>0.5;
        Xa(2,1:Kp) = Xa(1,1:Kp);
        
        % CRC Attachment
        Xa_crc = crcadd_m(Xa, poly);
        
        %SPC encoding: Compute the $X_{\hat{\mathcal{A}}}$
        ua = mod(Xa_crc * Gaa, 2);
        Xb = mod(ua * Gab, 2);
        
        %Set the bits respect to bit channel.
        codeword(info_index, :) = Xa_crc';
        codeword(frozen_index, :) = Xb';
        
        % bpsk modulation
        encode_sym = 1 - 2 * codeword;
        
        % add noise
        receive_sample = encode_sym + sigma * randn(size(encode_sym));
        
        %Compute the LLR; The LLR of the punctured bits should be set as 0.
        llr = 2*receive_sample/sigma^2;      

        %decoding of the block odd.
        [esti_u1, err1] = SPC_CASCL_decoder(llr(:,1), L, info_index, frozen_bits, Gaa, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        esti_Xa1 = mod(esti_u1' * Gaa, 2);
        decision_bits(1,:) = esti_Xa1;   %Store decision bits.
        
        %decoding of the block even.
        [esti_u2, err2] = SPC_CASCL_decoder(llr(:,2), L, info_index, frozen_bits, Gaa, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        esti_Xa2 = mod(esti_u2' * Gaa, 2);
        decision_bits(2,:) = esti_Xa2;   %Store decision bits.
        
        %If the CRC check is incorrect in the CA-SCL decoder, it should be
        %checked again.
        if err1
            err1 = sum(crccheck_m(esti_Xa1, poly));
        end
        
        if err2
            err2 = sum(crccheck_m(esti_Xa2, poly));
        end
        
        
        % situation 1: polar1 wrong, polar2 right;
        if err1 && ~err2
            
            ReSCL_oddWrong = ReSCL_oddWrong + 1;
            
            %Get extrinsic information from block even.
            llr(MUUB,1) = init_lr_max * (1 - 2 * decision_bits(2,1:Kp));
            
            %Second round of decoding.
            [esti_u1, ~] = SPC_CASCL_decoder(llr(:,1), L, info_index, frozen_bits, Gaa, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
            esti_Xa = mod(esti_u1' * Gaa, 2);
            decision_bits(1,:) = esti_Xa;
            
            %If correctly.
            if sum(crccheck_m(esti_Xa, poly))==0
                ReSCL_oddCorrect = ReSCL_oddCorrect + 1;
            end
        end
        
        %situation 2: polar1 right, polr2 wrong; In this case, no need for new round of decoding therein.
        if ~err1 && err2
            
            ReSCL_evenWrong = ReSCL_evenWrong + 1;
            
            %Get extrinsic information from block even.
            llr(MUUB,2) = init_lr_max * (1 - 2 * decision_bits(1,1:Kp));
            
            %Second round of decoding.
            [esti_u2, ~] = SPC_CASCL_decoder(llr(:,2), L, info_index, frozen_bits, Gaa, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
            esti_Xa2 = mod(esti_u2' * Gaa, 2);
            decision_bits(2,:) = esti_Xa2;
            
            %If correctly.
            if sum(crccheck_m(esti_Xa2, poly))==0
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
        % calculate BER and PER
        count1 = sum(decision_bits(1,:) ~= Xa_crc(1,:));
        if count1 ~= 0
            PerNum1 = PerNum1 + 1;
            BerNum1 = BerNum1 + count1;
        end
        count2 = sum(decision_bits(2,:) ~= Xa_crc(2,:));
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
% path = '../results/';
% filename = [path, 'PCM_SCL',num2str(L),'_N',num2str(N),'_R',num2str(R),'.mat'];
% save(filename)