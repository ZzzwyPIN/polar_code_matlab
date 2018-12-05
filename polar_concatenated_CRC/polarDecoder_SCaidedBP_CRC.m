clc
clear
% 基本参数设置
n = 8;  % 比特位数
R = 0.25;    % 码率
Ng = 8;
poly = [1 1 1 0 1 0 1 0 1];

SNR = 0:5;

init_lr_max = 3;    % limit the max LR of the channel to be with [-3 3]
max_iter = 40;
% block_num = 10000;

% 参数计算
snr = 10.^(SNR/10);
esn0 = snr * R;
init_max = init_lr_max * n;
if init_max > 30
    init_max = 30;
end
N = 2^n;
K = N*R;  % information bit length
Kp = N*R*0.25;  % Cascaded decoding length
k_f = N-K;% frozen_bits length
% source_block = 2*k-k1;
% frozen_block = 2*k_f;

% get information bits and concatenated bits
load('Pe_snr3p0db_2048_n_8.mat');   % load the channel information
[Ptmp, I] = sort(P);
info_index = sort(I(1:K));  % 挑选质量好的信道传输信息位
info_without_crc = sort(I(Ng+1:1:K));
frozen_index = sort(I(K+1:end));   % 传输冻结位的信道
% inter_index = sort(I(Ng+1:Ng+Kp)); % 相互信息位
inter_index = sort(I(K-Kp+1:K));
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
    ReBP_counter = 0;
    ReBP_correct = 0;
    iter = 0;
    while true
        iter = iter + 1;
        fprintf('\nNow iter: %2d\tNow SNR: %d', iter, SNR(i));
        source_bit1 = randi([0 1],1,K-Ng);
        source_bit2 = randi([0 1],1,K-Kp-Ng);
        [~,temp_index] = ismember(inter_index,info_without_crc);
        source_bit2 = insert_bit(source_bit1,source_bit2,temp_index,temp_index);
        source_crc_bit1 = crcadd(source_bit1,poly);
        source_crc_bit2 = crcadd(source_bit2,poly);
        encode_temp1 = rem(source_crc_bit1*Gi + frozen_bits*Gf,2);
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
            [lr_u,lr_x] = getBP_Parameter(receive_sample1,frozen_bits,frozen_index,n,init_max,sigma);
            for m = 1:length(temp_index)
                if decision_bits2(temp_index(m)) == 0
                    lr_u(reverse_index(n,info_index(temp_index(m)))) = init_max;
                else
                    lr_u(reverse_index(n,info_index(temp_index(m)))) = -init_max;
                end
            end
            decision_bits1 = polarBP_decoder(n,lr_u,lr_x,max_iter,info_index);
            ReBP_counter = ReBP_counter + 1;
            if sum(crccheck(decision_bits1,poly)) == 0
               ReBP_correct =  ReBP_correct + 1;
            end
        end
        
        % situation 2: polar1 right, polr2 wrong;
        if isempty(find(receive_crc_bits1,1)) && ~isempty(find(receive_crc_bits2,1))
            % get init LLR
            [lr_u,lr_x] = getBP_Parameter(receive_sample2,frozen_bits,frozen_index,n,init_max,sigma);
            for m = 1:length(temp_index)
                if decision_bits1(temp_index(m)) == 0
                    lr_u(reverse_index(n,info_index(temp_index(m)))) = init_max;
                else
                    lr_u(reverse_index(n,info_index(temp_index(m)))) = -init_max;
                end
            end
            decision_bits2 = polarBP_decoder(n,lr_u,lr_x,max_iter,info_index);
            ReBP_counter = ReBP_counter + 1;
            if sum(crccheck(decision_bits2,poly)) == 0
               ReBP_correct =  ReBP_correct + 1;
            end
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
        if (PerNum1 >= 100 && PerNum2 >= 100 && iter >= 10000)
            break;
        end
    end
    per(i) = (PerNum1+PerNum2)/(2*iter);
    ber(i) = (BerNum1+BerNum2)/(2*K-Kp)/iter;
    rs_coun(i) = ReBP_counter;
    rs_corr(i) = ReBP_correct;
end
fprintf('\nNow disp the Ber and Per');
fprintf('\nPer\t\tBer\t\tPer1\tBer1\tPer2\tBer2\tEbN0');
for i = 1:length(SNR)
    fprintf('\n%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%d',per(i),ber(i),per1(i),ber1(i),per2(i),ber2(i),SNR(i));
end
semilogy(SNR,per1,'b-*',SNR,ber1,'b-+',SNR,per2,'k-*',SNR,ber2,'k-+',SNR,per,'r-*',SNR,ber,'r-+');
xlabel('SNR in dB');
ylabel('BER and PER in dB');
title('Cascaded Polar Decoding');
legend('PER1','BER1','PER2','BER2','PER','BER');

% 
% 

%%%0.25 total ber
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR,berSC_ReSC_w,'r-d',SNR,berSC_ReBP,'m-p',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(0:4,berSC_ReSC_b,'r-s',0:4,berBP_ReSC,'k-h',0:4,berBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([0 4 6.1e-06 0.4])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')

%%%0.25 total per
semilogy(SNR,perSC,'b-*',SNR,perBP,'k-+',SNR,perSC_ReSC_w,'r-d',SNR,perSC_ReBP,'m-p',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
semilogy(0:4,perSC_ReSC_b,'r-s',0:4,perBP_ReSC,'k-h',0:4,perBP_ReBP,'b-p');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([0 4 2.2e-05 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)','BP-SC','BP-BP')

%%%0.25 compare
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR,ber,'r-d',0:4,berSC_SC_best,'r-s');
hold on
semilogy(SNR,perSC,'b--*',SNR,perBP,'k--+',SNR,per,'r--d',0:4,perSC_SC_best,'r--s');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([0 4 6.1e-06 0.4])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-SC(K_{p}best)')

%%%%0.25ber
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR,ber,'r-d',SNR,berCA_BP,'m-p',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(0:4,berSC_SC_best,'r--s');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([0 4 6.1e-06 0.4])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)')

%%%0.5ber
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR,berCA_SC,'r-d',SNR,berCA_BP,'m-p',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(0:4,berSC_SC_best,'r--s');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([0 4 9.4e-06 1])
grid on
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)')

%%%0.75ber
semilogy(SNR,berSC,'b-*',SNR,berBP,'k-+',SNR,ber,'r-d',SNR,berCA_BP,'m-p',SNR,berSCL8,'c-o',SNR,berSCL32,'b-^');
hold on
semilogy(0:4,berSC_SC_best,'r--s');
xlabel('SNR in dB');
ylabel('Bit Error rate');
axis([1 5 1.0e-05 1])
legend('SC','BP','SC-SC(K_{p}worst)','SC-BP(K_{p}worst)','SCL8','SCL32','SC-SC(K_{p}best)')
grid on

%%%per0.25
semilogy(SNR,perSC,'b-*',SNR,perBP,'k-+',SNR,perCA_SC,'r-d',SNR,perCA_BP,'m-p',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
xlabel('SNR in dB');
ylabel('BER and PER in dB');
axis([0 4 2.2e-05 1])
grid on
legend('SC','BP','CA-SC','CA-BP','SCL8','SCL32')

%%%per0.5
semilogy(SNR,perSC,'b-*',SNR,perBP,'k-+',SNR,perCA_SC,'r-d',SNR,perCA_BP,'m-p',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
xlabel('SNR in dB');
ylabel('BER and PER in dB');
axis([0 4 4.4e-05 1])
grid on
legend('SC','BP','CA-SC','CA-BP','SCL8','SCL32')

%%%per0.75
semilogy(SNR,perSC,'b-*',SNR,perBP,'k-+',SNR,perCA_SC,'r-d',SNR,perCA_BP,'m-p',SNR,perSCL8,'c-o',SNR,perSCL32,'b-^');
hold on
xlabel('SNR in dB');
ylabel('BER and PER in dB');
axis([1 4 6.00e-05 1])
grid on
legend('SC','BP','CA-SC','CA-BP','SCL8','SCL32')