clc
clear

% 基本参数设置
n = 8;  % 比特位数
R = 0.5;    % 码率

SNR = -1:5;

init_lr_max = 3;    % limit the max LR of the channel to be with [-3 3]
max_iter = 30;
block_num = 10;

% 参数计算
snr = 10.^(SNR/10);
init_max = init_lr_max * n;
if init_max > 30
    init_max = 30;
end
N = 2^n;
K = floor(N*R);  % information bit length
k_f = N - K;

% get information bits and concatenated bits
load('Pe_snr3p0db_2048_n_8.mat');   % load the channel information
[Ptmp, I] = sort(P);
Info_index = sort(I(K:-1:1));  % 挑选质量好的信道传输信息位
Frozen_index = sort(I(end:-1:K+1));   % 传输冻结位的信道

% get generate matrix
G = encoding_matrix(n);
Gi = G(Info_index,:);
Gf = G(Frozen_index,:);
frozen_bits = zeros(1,k_f);
rng('shuffle')
for i = 1:length(SNR)
    sigma = (1/snr(i))^0.5;
    % set PER and BER counter
    PerNumSC = 0;
    BerNumSC = 0;
    PerNumBP = 0;
    BerNumBP = 0;
    for iter = 1:block_num
        fprintf('\nNow iter: %2d\tNow SNR: %d', iter, SNR(i));
        source_bit = randi([0 1],1,K);
        encode_temp = rem(source_bit*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp = (-1).^(encode_temp + 1);
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        % SC decoder
        receive_bits_SC = polarSC_decoder(n,receive_sample,snr(i),Frozen_index,frozen_bits,Info_index);
        
        % BP decoder
        lr_x = -2*receive_sample./(sigma^2);
        % decoding follow
        lr_u = zeros(1,N); % save send sample LR in each iteration
        lr_u(reverse_index(n,Frozen_index)) = init_max;
        
        receive_bits_BP = polarBP_decoder(n,lr_u,lr_x,max_iter,Info_index);
        
        % calculate BER and PER
        countSC = sum(receive_bits_SC ~= source_bit);
        if countSC ~= 0
            PerNumSC = PerNumSC + 1;
            BerNumSC = BerNumSC + countSC;
        end 
        
        countBP = sum(receive_bits_BP ~= source_bit);
        if countSC ~= 0
            PerNumBP = PerNumBP + 1;
            BerNumBP = BerNumBP + countBP;
        end
    end
    perSC(i) = PerNumSC/block_num;
    berSC(i) = BerNumSC/(K*block_num);
    perBP(i) = PerNumBP/block_num;
    berBP(i) = BerNumBP/(K*block_num);
end
semilogy(SNR,perSC,'k-+',SNR,berSC,'k-*',SNR,perBP,'b-+',SNR,berSC,'b-*')
xlabel('SNR in dB');
ylabel('BER and PER in dB');
title('Cascaded Polar Decoding');
legend('perSC','berSC','perBP','berSC');