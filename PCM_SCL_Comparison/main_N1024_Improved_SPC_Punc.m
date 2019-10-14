%Improved PCM using SPC ecoding and decoding. And the punctured is used to
%keep the code rate the same as that of the stand-alone SC.

clc
clear

% 基本参数设置
n = 10;  % 比特位数
Ng = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
L = 8;   %SCL List
K = 528; %the number of information bits of the underlying blocks
Kp = 40; %the number of mutual bits
block_num = 2;
init_max = 3;
SNR = 0.5:0.5:3.5;

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
    PerNum = 0;
    BerNum = 0;
    iter = 0;
    %counter the number of Re-SC decoding
    % 以下参数用来记录每个SNR点，论文中提到的case1-case4发生次数
    ReSCL_oddWrong = 0;
    ReSCL_oddCorrect = 0;    %odd block have correct new rounds of SCL decoding
    oddCorrect = 0; %Counter the number of the block odd decode correctly in the first round of decoding.
    evenCorrect = 0;%Counter the number of the block even decode correctly in the first round of decoding, when block odd is right.
    AllRight = 0;
    AllWrong = 0;
    
    while true 
        
        iter = iter + 1;
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum: %2d\tNow Error Bits: %2d', iter, SNR(i), PerNum, BerNum);
        
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
        
        %Puncture: Set the punctured symbols as 0.
        encode_sym(MUUB, 2) = 0;
        
        % add noise
        receive_sample = encode_sym + sigma * randn(size(encode_sym)).*abs(encode_sym);
        
        %Compute the LLR; The LLR of the punctured bits should be set as 0.
        llr = 2*receive_sample/sigma^2;      

        %decoding of the block odd.
        [esti_u, err1] = SPC_CASCL_decoder(llr(:,1), L, info_index, frozen_bits, Gaa, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        esti_Xa = mod(esti_u' * Gaa, 2);
        
        %If the CRC check is incorrect in the CA-SCL decoder, it should be
        %checked again.
        if err1
            err1 = sum(crccheck_m(esti_Xa, poly));
        end
        
        decision_bits(1,:) = esti_Xa;   %Store decision bits.
        
        %If the block one decode correctly, it can provide block 2
        %correct mutual information bits. We can regard those bits as
        %extrinsic information, so their LLR can be set as inf * (1 - 2 * u_i).
        if ~err1
            llr(MUUB, 2) = init_lr_max * (1 - 2 * decision_bits(1,1:Kp));
            oddCorrect = oddCorrect + 1;
        else
            llr(MUUB, 2) = llr(MUUB, 1);
        end
        
        [esti_u, err2] = SPC_CASCL_decoder(llr(:,2), L, info_index, frozen_bits, Gaa, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
        esti_Xa = mod(esti_u' * Gaa, 2);
        if err2
            err2 = sum(crccheck_m(esti_Xa,poly));
        end
        decision_bits(2,:) = esti_Xa;   %Store decision bits
        
        % situation 1: polar1 wrong, polar2 right;
        if err1 && ~err2
            
            ReSCL_oddWrong = ReSCL_oddWrong + 1;
            
            %Get extrinsic information from block even.
            llr(MUUB,1) = init_lr_max * (1 - 2 * decision_bits(2,1:Kp));
            
            %Second round of decoding.
            [esti_u, ~] = SPC_CASCL_decoder(llr(:,1), L, info_index, frozen_bits, Gaa, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
            esti_Xa = mod(esti_u' * Gaa, 2);
            decision_bits(1,:) = esti_Xa;
            
            %If correctly.
            if sum(crccheck_m(esti_Xa, poly))==0
                ReSCL_oddCorrect = ReSCL_oddCorrect + 1;
            end
        end
        
        %situation 2: polar1 right, polr2 wrong; In this case, no need for new round of decoding therein.
        if ~err1 && err2
            evenCorrect = evenCorrect + 1;
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
        count = sum(decision_bits ~= Xa_crc, 'all');
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end
        
        if (PerNum>=100 && iter>=10000)
            break;
        end
        
        if iter >= 10000000
           break; 
        end
        
    end
    iterNum(i) = iter;
    per(i) = PerNum/iter;
    ber(i) = BerNum/2/K/iter;
    rs_oddwrong(i) = ReSCL_oddWrong;
    rs_oddcorr(i) = ReSCL_oddCorrect;
    evencorr(i) = evenCorrect;
    oddcorr(i) = oddCorrect;
    all_right(i) = AllRight;
    all_wrong(i) = AllWrong;
end

% record the results
% path = '../results/';
% filename = [path, 'PCM_SCL',num2str(L),'_N',num2str(N),'_R',num2str(R),'.mat'];
% save(filename)