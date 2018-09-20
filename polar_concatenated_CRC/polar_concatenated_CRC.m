clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.5;    % 码率
Ng = 8;
poly = [1 1 1 0 1 0 1 0 1];

SNR = -2:5;

init_lr_max = 3;    % limit the max LR of the channel to be with [-3 3]
max_iter = 30;
block_num = 10000;

% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
init_max = init_lr_max * n;
if init_max > 30
    init_max = 30;
end
N = 2^n;
K = N*R;  % information bit length
k = N*R*R*R;  % Cascaded decoding length
k_f = N-K;% frozen_bits length
% source_block = 2*k-k1;
% frozen_block = 2*k_f;

% get information bits and concatenated bits
load('Pe_snr3p0db_2048_n_8.mat');   % load the channel information
[Ptmp, I] = sort(P);
info_index = sort(I(K:-1:1));  % 挑选质量好的信道传输信息位
frozen_index = sort(I(end:-1:K+1));   % 传输冻结位的信道
bad_info_index = sort(I(K:-1:K-k+1)); % 级联解码信息位
% best_info_index = sort(I(k:-1:1));

% get generate matrix
G = encoding_matrix(n);
Gi = G(info_index,:);
Gf = G(frozen_index,:);
frozen_bits = randi([0 1],1,k_f);
rng('shuffle');
for i = 1:length(SNR)
    sigma = (2*esn0(i))^(-0.5);
    % set PER and BER counter
    PerNum1 = 0;
    BerNum1 = 0;
    PerNum2 = 0;
    BerNum2 = 0;
    for iter = 1:block_num
        fprintf('\nNow iter: %2d\tNow SNR: %d', iter, SNR(i));
        source_bit1 = randi([0 1],1,K-Ng);
        source_bit2 = randi([0 1],1,K-k-Ng);
        source_crc_bit1 = crcadd(source_bit1,poly);
        encode_temp1 = rem(source_crc_bit1*Gi + frozen_bits*Gf,2);
        [~,temp_index] = ismember(bad_info_index,info_index);
        source_bit2 = insert_bit(source_crc_bit1,source_bit2,temp_index,temp_index);
        source_crc_bit2 = crcadd(source_bit2,poly);
        encode_temp2 = rem(source_crc_bit2*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp1 = (-1).^(encode_temp1 + 1);
        encode_temp2 = (-1).^(encode_temp2 + 1);
        % add noise
        receive_sample1 = encode_temp1 + sigma * randn(size(encode_temp1));
        receive_sample2 = encode_temp2 + sigma * randn(size(encode_temp2));
        
        % SC decoder follow
        decision_bits1 = polarSC_decoder(n,receive_sample1,sigma,frozen_index,frozen_bits,info_index);
        decision_bits2 = polarSC_decoder(n,receive_sample2,sigma,frozen_index,frozen_bits,info_index);
        
        % CRC check follow
        receive_crc_bits1 = crccheck(decision_bits1,poly);
        receive_crc_bits2 = crccheck(decision_bits2,poly);
        % crc Check Result：If only one polar is uncorrect,then using BP
        % decoder with some concatenated bits extrasinc information.
        
        % situation 1: polar1 wrong, polr2 right;
        if ~isempty(find(receive_crc_bits1,1)) && isempty(find(receive_crc_bits2,1))
            % get init LLR
            lr_x = -2*receive_sample1./(sigma^2);
            % decoding follow
            lr_u = zeros(1,N); % save send sample LR in each iteration
            frozen_index_0 = frozen_bits == 0;
            frozen_index_1 = frozen_bits == 1;
            lr_u(reverse_index(n,frozen_index(frozen_index_0))) = init_max;
            lr_u(reverse_index(n,frozen_index(frozen_index_1))) = -init_max;
            for m = 1:length(temp_index)
                if decision_bits2(temp_index(m)) == 0
                    lr_u(reverse_index(n,info_index(temp_index(m)))) = init_max;
                else
                    lr_u(reverse_index(n,info_index(temp_index(m)))) = -init_max;
                end
            end
            decision_bits1 = polarBP_decoder(n,lr_u,lr_x,max_iter,info_index);
        end
        
        % situation 2: polar1 right, polr2 wrong;
        if isempty(find(receive_crc_bits1,1)) && ~isempty(find(receive_crc_bits2,1))
            % get init LLR
            lr_x = -2*receive_sample2./(sigma^2);
            % decoding follow
            lr_u = zeros(1,N); % save send sample LR in each iteration
            frozen_index_0 = frozen_bits == 0;
            frozen_index_1 = frozen_bits == 1;
            lr_u(reverse_index(n,frozen_index(frozen_index_0))) = init_max;
            lr_u(reverse_index(n,frozen_index(frozen_index_1))) = -init_max;
            for m = 1:length(temp_index)
                if decision_bits1(temp_index(m)) == 0
                    lr_u(reverse_index(n,info_index(temp_index(m)))) = init_max;
                else
                    lr_u(reverse_index(n,info_index(temp_index(m)))) = -init_max;
                end
            end
            decision_bits2 = polarBP_decoder(n,lr_u,lr_x,max_iter,info_index);
        end
        
        % situation 3 and 4: polar1 and polar2 are both right or wrong
        % we have no salution.
        
        % calculate BER and PER
        count1 = sum(decision_bits1 ~= source_crc_bit1);
        if count1 ~= 0
            PerNum1 = PerNum1 + 1;
            BerNum1 = BerNum1 + count1;
        end
        count2 = sum(decision_bits2 ~= source_crc_bit2);
        if count2 ~= 0
            PerNum2 = PerNum2 + 1;
            BerNum2 = BerNum2 + count2;
        end
    end
    per1(i) = PerNum1/block_num;
    per2(i) = PerNum2/block_num;
    ber1(i) = BerNum1/(K*block_num);
    ber2(i) = BerNum2/(K*block_num);
    per(i) = (per1(i)+per2(i))/2;
    ber(i) = (BerNum1+BerNum2)/(2*K-k)/block_num;
end
fprintf('\nNow disp the Ber and Per');
fprintf('\nPer\t\tBer\t\tPer1\tBer1\tPer2\tBer2\tEbN0');
for i = 1:length(SNR)
    fprintf('\n%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%d',per(i),ber(i),per1(i),ber1(i),per2(i),ber2(i),SNR(i));
end
semilogy(SNR,per1,'b-*',SNR,ber1,'b-+',SNR,per2,'k-*',SNR,ber2,'k-+',SNR,per,'r-*',SNR,ber,'r-+');
hold on
semilogy(SNR,perSC,'g-*',SNR,berSC,'g-+',SNR,perBP,'k-^',SNR,berBP,'k-o');
xlabel('SNR in dB');
ylabel('BER and PER in dB');
title('Cascaded Polar Decoding');
legend('PER1','BER1','PER2','BER2','PER','BER','perSC','berSC','perBP','berBP');