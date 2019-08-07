clc
clear

% 基本参数设置
n = 10;  % 比特位数 
N = 2^n;
Ng = 12;
poly = [1 1 1 1 1 0 0 0 1 0 0 1 1];


K = 512+12;
Kinfo = 512;
Kp = 36;

% Kpure = K-Ng-Kp;
k_f = N-K;% frozen_bits length
% R = K/N;

SNR = 3;
% 参数计算
snr = 10.^(SNR/10);
% esn0 = snr * R;
% sigma = (2*esn0)^(-0.5);

lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);




rng('shuffle');
for i = 1:length(Kp)
    
    Kpure = K-Ng-Kp(i);
    R = (K-Ng-Kp(i)/2)/N;
    
    load('Pe_N1024_snr3_R5.mat');
    [~, I] = sort(P);
    pure_index = I(1:Kpure);
    inter_index = I(K-Kp(i)+1:K);
    crc_index = I(Kpure+1:K-Kp(i));
    info_index = [pure_index inter_index crc_index];
    frozen_index = I(K+1:end);   % 传输冻结位的信道
    
    for ii = 1:length(SNR)
        
        
        esn0 = snr(ii)*R;
        sigma = (2*esn0)^(-0.5);
        
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
        ReSC_oddCorrect = 0;    %odd block have correct new rounds of SCL decoding
        ReSC_evenCorrect = 0;   %even block have correct new rounds of SCL decoding
        AllRight = 0;
        AllWrong = 0;

        %以下参数用来记录miss check of CRC

        missCRCodd = 0;
        missCRCeven = 0;
        missAll = 0;

        while true 


            flag1 = 0;
            flag2 = 0;
            flag3 = 0;

            iter = iter + 1;
            % reset the frozen bits and mutual bits
            frozen_bits = ones(N,1);
            mutual_bits = zeros(N,1);
            frozen_bits(info_index) = 0;

            if mod(iter,1000)==0
                fprintf('\nNow iter: %2d\tNow Kp: %d\tNow SNR: %2d\tNow PerNum1: %2d\tNow PerNum2: %2d\tNow Error Bits: %2d',iter,Kp(i),SNR(ii),PerNum1,PerNum2,BerNum1+BerNum2);
            end

            source_bit1 = rand(1,Kinfo)>0.5;
            source_bit2 = zeros(1,Kinfo);
            source_bit2(1:Kpure) = rand(1,Kpure)>0.5;
    %         [~,temp_index] = ismember(inter_index,info_without_crc);
    %         source_bit2 = insert_bit(source_bit1,source_bit2,temp_index,temp_index);
            source_bit2(Kpure+1:Kinfo) = source_bit1(Kpure+1:Kinfo);
            source_crc_bit1 = crcadd(source_bit1,poly);
            source_crc_bit2 = crcadd(source_bit2,poly);

            u1 = zeros(N, 1);
            u2 = zeros(N, 1);
            u1(info_index) = source_crc_bit1;
            u2(info_index) = source_crc_bit2;
    %         u1 = getU(u1,pure_index,inter_index,crc_index,source_crc_bit1);
    %         u2 = getU(u2,pure_index,inter_index,crc_index,source_crc_bit2);

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


            decision_bits1 = SC_decoder(llr1, info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
            decision_bits2 = SC_decoder(llr2, info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);

    %         [decision_bits1, err1] = CASCL_decoder(llr1, L, K, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);
    %         [decision_bits2, err2] = CASCL_decoder(llr2, L, K, frozen_bits, poly, lambda_offset, llr_layer_vec, bit_layer_vec);

            err1 = sum(crccheck(decision_bits1',poly))==0;
            err2 = sum(crccheck(decision_bits2',poly))==0;
            % crc Check Result：If only one polar is uncorrect,then using BP
            % decoder with some concatenated bits extrasinc information.

            % situation 1: polar1 wrong, polr2 right;
            if ~err1 && err2
                flag2 = 1;
                % modify polar1 frozen_index frozen_bits info_index
                ReSC_oddWrong = ReSC_oddWrong + 1;
                % modify polar1 frozen_index frozen_bits info_index
                frozen_bits(inter_index) = 2;
                mutual_bits(inter_index) = decision_bits2(Kpure+1:Kinfo);
                decision_bits1 = SC_decoder(llr1, info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                if sum(decision_bits1(1:Kinfo)' ~= source_bit1) == 0
                    ReSC_oddCorrect = ReSC_oddCorrect + 1;
                end
            end

            % situation 2: polar1 right, polr2 wrong;
            if err1 && ~err2
                flag1 = 1;
                ReSC_evenWrong = ReSC_evenWrong + 1;
                % modify polar1 frozen_index frozen_bits info_index
                frozen_bits(inter_index) = 2;
                mutual_bits(inter_index) = decision_bits1(Kpure+1:Kinfo);
                decision_bits2 = SC_decoder(llr2, info_index, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, mutual_bits);
                if sum(decision_bits2(1:Kinfo)' ~= source_bit2) == 0
                    ReSC_evenCorrect = ReSC_evenCorrect + 1;
                end
            end

            if ~err1 && ~err2
                AllWrong = AllWrong + 1;
            end

            if err1 && err2
                flag3 = 1;
                AllRight = AllRight + 1;
            end


            % situation 3 and 4: polar1 and polar2 are both right or wrong
            % we have no salution.

            % calculate BER and PER
            count1 = sum(decision_bits1(1:Kinfo)' ~= source_bit1);
            if count1 ~= 0
                if flag1
                    missCRCodd = missCRCodd +1;
                end
                PerNum1 = PerNum1 + 1;
                BerNum1 = BerNum1 + count1;
            end
            count2 = sum(decision_bits2(1:Kinfo)' ~= source_bit2);
            if count2 ~= 0
                if flag2
                    missCRCeven = missCRCeven +1;
                end
                PerNum2 = PerNum2 + 1;
                BerNum2 = BerNum2 + count2;
            end

            if ((count1 ~= 0 || count2 ~=0)&& flag3)
               missAll = missAll +1; 
            end

            if (PerNum1>=100 && PerNum2>=100 && iter>=10000)
                break;
            end
            if iter>=10000000
               break; 
            end

        end
        enditerNum(i,ii) = iter;
        per(i,ii) = (PerNum1+PerNum2)/(2*iter);
        ber(i,ii) = (BerNum1+BerNum2)/(2*Kinfo-Kp(i))/iter;
        rs_oddwrong(i,ii) = ReSC_oddWrong;
        rs_evenwrong(i,ii) = ReSC_evenWrong;
        rs_oddcorr(i,ii) = ReSC_oddCorrect;
        rs_evencorr(i,ii) = ReSC_evenCorrect;
        all_right(i,ii) = AllRight;
        all_wrong(i,ii) = AllWrong;

        missodd(i,ii) = missCRCodd;
        misseven(i,ii) = missCRCeven;
        missall(i,ii) = missAll;
    end
end
success_rate = (rs_oddcorr+rs_evencorr)./(rs_oddwrong+rs_evenwrong);
% record simulation results
% path = './results/';
% filename = [path, 'PCM_FastSC_N',num2str(N),'_R',num2str(R),'.mat'];
% save(filename)