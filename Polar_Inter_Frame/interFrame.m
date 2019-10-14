clc
clear

% 基本参数设置
n = 10;  % 比特位数 
N = 2^n;
Ng = 16;
L = 16;
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
m = 40; % The number of the MUUB and MRFB.
M = 6; %The number of encoding block.

SNR = 1.5:0.5:4;

K = 528;
Kpure = K - Ng;
R = Kpure/N;

snr = 10.^(SNR/10);
esn0 = snr * R;

load('Pe_N1024_snr2.mat');
[~, I] = sort(P);
crc_index = I(1:Ng);
pure_info_index = I(Ng+1:K-m); %Pure information index.
MUUB = I(K-m+1:K); %Most unreliable unfrozen bits.
MRFB = I(K+1:K+m); %Most reliable frozen bits.
frozen_index = I(K+m+1:end);   %Pure frozen index.
info_index = [MRFB MUUB pure_info_index crc_index];

lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

rng('shuffle');
for i = 1:length(SNR)
    
    sigma = (2*esn0(i))^(-0.5);
    iter = 0;
    BerNum = 0;
    PerNum = 0;
    while true
        
        iter = iter + 1;
        
        frozen_bits = ones(N,1);
        frozen_bits([MUUB pure_info_index crc_index]) = 0;
        mutual_bits = zeros(N,1);
        u = zeros(N,M);
        decision_bits = zeros(K+m,M);
        
        % Set the indicator g
        g = 0;
        
        fprintf('\nNow iter: %2d\tNow SNR: %2d\tNow PerNum: %2d\tNow BerNum: %2d',iter, SNR(i), PerNum, BerNum); 
        
        % Encoding: u and encode_temp are the reverse bits vec and
        % codeword, respectively.
        source_bits = randn(M, K-Ng)>0.5;
        
        
%         source_bits = randn(1,K-Ng)>0.5;    %The source bits of the 1st CB.
%         source_bits_crc(1,:) = crcadd_m([zeros(1,m) source_bits], poly);
%         info_index = [MRFB MUUB pure_info_index];
%         u(info_index,1) = source_bits_crc(1,:);
%         
%         % encoder following for the i=2 to end.
%         for num = 2:M
%             source_bits = randn(1,K-Ng)>0.5;
%             source_bits_crc(num,:) = crcadd_m([source_bits_crc(num-1,m+1:2*m) source_bits], poly);
%             u(info_index,num) = source_bits_crc(num,:);
%         end
%         encode_temp = interFrameEncoder(u, lambda_offset, llr_layer_vec);
        
        %BPSK
        encode_temp = 1 - 2 * encode_temp;
        
        % add noise
        receive_sample = encode_temp + sigma * (randn(size(encode_temp)));
        
        % Decoding
        LLR = (2*receive_sample)/sigma^2;
        decision_bits(:,1) = CASCL_decoder(LLR(:,1), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
%         decision_bits(:,1) = SC_decoder(LLR(:,1), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        f = sum(crccheck_m(decision_bits(:,1)', poly));
        for n = 2:M
            if f
                frozen_bits(MRFB) = 0;
                decision_bits(:,n) = CASCL_decoder(LLR(:,n), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
                f = sum(crccheck(decision_bits(:,n)', poly));
                if ~f
                    if ~g && n>2
                        frozen_bits([MRFB MUUB]) = 2;
                        mutual_bits([MRFB MUUB]) = [source_bits_crc(n-2, m+1:2*m) source_bits_crc(n, 1:m)];
                        decision_bits(:,n-1) = CASCL_decoder(LLR(:,n-1), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec,mutual_bits);
%                         decision_bits(:,n-1) = SC_decoder(LLR(:,n-1), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                    elseif ~g && n==2
                        frozen_bits(MUUB) = 2;
                        mutual_bits(MUUB) = source_bits_crc(n, 1:m);
                        decision_bits(:,n-1) = CASCL_decoder(LLR(:,n-1), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
%                         decision_bits(:,n-1) = SC_decoder(LLR(:,n-1), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                    else
                        g = 0;
                    end
                else
                    g = 1;
                end
            else
                frozen_bits(MRFB) = 2;
                mutual_bits(MRFB) = source_bits_crc(n-1, m+1:2*m);
                decision_bits(:,n) = CASCL_decoder(LLR(:,n), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
%                 decision_bits(:,n) = SC_decoder(LLR(:,n), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                f = sum(crccheck(decision_bits(:,n)', poly));
            end
            
            %reset the frozen_bits
            frozen_bits(:) = 1;
            frozen_bits([MUUB pure_info_index]) = 0;
        end
        
        %Counter the Error bits and blocks
        counter = sum(source_bits_crc ~= decision_bits', 'all');
        if counter ~= 0
            PerNum = PerNum + 1;
        end
        BerNum = BerNum + counter;
        
        if (PerNum >= 100 && iter >= 10000)
            break;
        end
        
        if (iter >= 10000000)
            break;
        end
    end
    ber(i) = BerNum/(K*M*iter);
    per(i) = PerNum/iter;
end