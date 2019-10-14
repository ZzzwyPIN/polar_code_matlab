clc
clear

% 基本参数设置
n = 10;  % 比特位数
SNR =  2;
Ng = 16;    %CRC bits number.
poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1]; %CRC generation polynomial equation.
L = 8;  %CA_SCL decoder with L lists
Lw = 10;
Rc = 0.15;
K = 558;
init_max = 3;

% 参数计算
N = 2^n;
init_lr_max = init_max * n;
Kc = floor(K*Rc);
Ku = K - 2*Kc;
R = (Lw*K-(Lw+1)*Kc)/(Lw*N-(Lw+1)*Kc);
R = floor(R*10)/10;
snr = 10.^(SNR/10);
esn0 = snr * R;
k_f = N-K;% frozen_bits length

% bit channel sort
load('Pe_N1024_snr2.mat');
[~, I] = sort(P);
crc_index = I(1:Ng);
pure_info_index = I(Ng+1:K-2*Kc);
lenp = length(pure_info_index);
MUUB = I(K-2*Kc+1:K);
Bh = MUUB(1:2:end);
Bt = MUUB(2:2:end);
frozen_index = I(K+1:end);

info_index = [Bh pure_info_index Bt crc_index];
% flag indicate if the bit channel carrys a info or frozen bit.
frozen_bits = ones(N,1);
frozen_bits(info_index) = 0;% 挑选质量好的信道传输信息位

%Generation matrix
G = spc_encoding(n);
Gaa = G(info_index, info_index);
Gab = G(info_index, frozen_index);

% Some parameters in the decoding process.
lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

rng('shuffle');
for i = 1:length(SNR)
    
    %Compute the noise varience
    sigma = (2*esn0(i))^(-0.5);
    
    % set PER and BER counter
    PerNum = 0;
    iter = 0;
    last = 0;
    back_err = 0;
    
    while true 
        iter = iter + 1;
        fprintf('\nNow iter: %2d\tNow SNR: %d\tNow PerNum: %2d',iter,SNR(i),PerNum);
        %Arrangement of the source bits
        codeword = zeros(N, Lw);
        Xa = rand(Lw,K-Ng)>0.5;
        Xa(1,1:Kc) = 0;
        Xa(2:end,1:Kc) = Xa(1:end-1,K-Ng-Kc+1:end);
        Xa(Lw,K-Ng-Kc+1:end) = 0;
        Xa_crc = crcadd_m(Xa, poly);
        
        %Compute the u
        ua = mod(Xa_crc * Gaa, 2);
        
        %Compute the codeword of the frozen bits
        Xb = mod(ua * Gab, 2);
        codeword(info_index,:) = Xa_crc';
        codeword(frozen_index,:) = Xb';
        
        % bpsk modulation
        encode_sym = 1 - 2 * codeword;
        
        %Puncture: Set the punctured symbols as 0.
        encode_sym(Bh,:) = 0;
        encode_sym(Bt,Lw) = 0;
        
        % add noise
        receive_sample = encode_sym + sigma * randn(size(encode_sym)).*abs(encode_sym);
        
        %Compute the LLR; The LLR of the punctured bits should be set as 0.
        llr = 2*receive_sample/sigma^2;
        
        %punctured the tail bits of the last block and the head bits of the
        %first block.
        llr(Bh,1) = init_lr_max;
        llr(Bt,Lw) = init_lr_max;
        
        %Decoding:look-back and go-back
        %Initialize:
        d = zeros(1,Lw);
        decision_bits = zeros(length(info_index), Lw);
        
        % decoding of the 1st CB. Feed-forward.
        [esti_bits, ~] = CASCL_decoder(llr(:,1), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
%         esti_bits = SC_decoder(llr(:,1), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        esti_u = mod(esti_bits'* Gaa, 2);
        decision_bits(:,1) = esti_u;
        d(1) = sum(crccheck_m(esti_u, poly))==0;   %set d = 1 if the esti bits pass the crc check.
        m = 2;  %the loop counter
        M = 0;  %the break point of the feed-forward decoding
        ff_counter = 0; %Feed-forward incorrect decoding
        while (m <= Lw)
            
            mutual_bits = zeros(N,1);
            % flag indicate if the bit channel carrys a info or frozen bit.
            frozen_bits = ones(N,1);
            frozen_bits(info_index) = 0;% 挑选质量好的信道传输信息位

            if d(m) && ~d(m-1)  % feed-back decoding
                M = m;
                m = m - 1;
                
                %Compute the llr of the m th CB.
                %The LLR of the index Bh has been computed in the
                %feed-forward decoding.
                llr(Bt,m) = init_lr_max * (1 - 2*decision_bits(1:Kc,m+1));
                mutual_bits(Bt) = decision_bits(1:Kc,m+1);
                frozen_bits(Bt) = 2;
                [esti_bits, ~] = CASCL_decoder(llr(:,m), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
%                 esti_bits = SC_decoder(llr(:,m), info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
                esti_u = mod(esti_bits'* Gaa, 2);
                d(m) = sum(crccheck_m(esti_u, poly))==0;   %set d = 1 if the esti bits pass the crc check.
                decision_bits(:,m) = esti_u;
                             
                if m == 1
                    if ~d(m)
                        back_err = back_err + 1;
                        break;
                    else
                        m = M + 1;
                        ff_counter = 0;
                    end
                else
                    if ~d(m)
                        back_err = back_err + 1;
                        break;
                    elseif d(m) && d(m-1)
                        m = M + 1;
                        ff_counter = 0;
                    end
                end
                
                
            else % feed-forward decoding
                if d(m-1)
                    llr(Bh,m) = init_lr_max * (1 - 2*decision_bits(Kc+lenp+1:2*Kc+lenp,m-1));
                    X_temp = decision_bits(Kc+lenp+1:2*Kc+lenp,m-1);
                    
                    frozen_bits(Bh) = 2;
                    [esti_bits, ~] = CASCL_decoder(llr(:,m), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                else
                    llr(Bh,m) = llr(Bt,m-1);
                    [esti_bits, ~] = CASCL_decoder(llr(:,m), L, info_index, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
                end
                esti_u = mod(esti_bits'* Gaa, 2);
                d(m) = sum(crccheck_m(esti_u, poly))==0;   %set d = 1 if the esti bits pass the crc check.
                decision_bits(:,m) = esti_u;
                
                if ~d(m)
                    ff_counter = ff_counter + 1;
                    if m==Lw
                        last = last + 1;
                    end
                end
                
                if ff_counter >= Lw
                    break;
                end
                
                if d(m) && ~d(m-1)
                    continue;
                else
                    m = m + 1;
                end
            end
        end
        
        count = sum(d);
        
        if count ~= Lw
            PerNum = PerNum + 1;
        end

        if (PerNum>=100 && iter>=10000)
            break;
        end
        
        if (iter >= 10000000)
           break; 
        end
        
        
    end    
    iterNum(i) = iter;
    per(i) = PerNum/iter;
end

% recording the results
% path = './results/';
% filename = [path, 'Polar_FastSC_N',num2str(N),'_R',num2str(R),'.mat'];
% save(filename)