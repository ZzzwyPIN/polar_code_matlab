% BP decoding of polar codes
% Dr Jan. 6, 2017
% Zzzwy modify in 2018/09/18
clc;
clear;

% diary on
% TOP_SCALE = 0;
n = 8;
N = 2^n;
R = 0.5; 
block_num = 10000;
SNR = -1:5;
snr = 10.^(SNR/10);
g_bp_init_max = 3;  % limit the max LR of the channel to be within [-3, 3]
init_max = g_bp_init_max * (n); % used to indicate max LR of frozen bits
%g_bp_inter_max = g_bp_max + n*g_bp_init_max; % intermediate max LR
if (init_max > 30)
   init_max = 30;
end
max_iter = 30;
K = floor(N*R); %infromation bits in one block
k_f = N - K;    %frozen_bits in one block
G = encoding_matrix(n); %产生G矩阵   F@(i)F     (1;n+1;2;n+2;3;n+3;......;n-1;2n-1;n;2n)


% source_bits=zeros(1,k);
frozen_bits = randi([0,1],1,N-K);
encoded_bits = zeros(1,N);  %极化码编码矩阵
received_bits = zeros(1,N);
load('Pe_snr0p0db_2048_n_8.mat');
% (2) find cut-off Z value below which the rows of G will be used for encoding
[Ptmp, I] = sort(P);
info_index = sort(I(K:-1:1));  % 挑选质量好的信道传输信息位
frozen_index = sort(I(end:-1:K+1));   % 传输冻结位的信道
Gi = G(info_index,:);  % 选取G矩阵coding行数
% (6) form the encoding matrix for frozen bits
Gf = G(frozen_index,:);  % 选取G矩阵frozen行数
rng('shuffle');
ebn0 = K/N*snr;
for i_aw = 1:length(SNR)
    sigma = (2*ebn0(i_aw))^(-0.5);
    ber = 0;
    per = 0;
    % BPSK modulation
    for r_idx = 1:block_num
        fprintf('\nNow iter: %2d\tNow SNR: %d',r_idx, SNR(i_aw));
        noise = sigma * randn(1,N);
        source_bits = randi([0,1],1,K);
        encoded_bits = source_bits * Gi + frozen_bits * Gf;
        idx_0 = find(mod(encoded_bits,2) ==  0);
        idx_1 = find(mod(encoded_bits,2) ==  1);
        encoded_bits(idx_0) = -1*ones(1,length(idx_0));
        encoded_bits(idx_1) = ones(1,length(idx_1));
        received_sample = encoded_bits + noise;%noise is added
        decision_bits = polar_bp_decoderM(received_sample,sigma,info_index,frozen_bits,frozen_index,max_iter,init_max,n);
        
        count = sum(decision_bits ~= source_bits);
        if count ~= 0
           ber = ber + count;
           per = per +1;
        end
    end
    BER_BP(i_aw)= ber/(K * block_num);
    PER_BP(i_aw)= per / block_num;
end
semilogy(SNR,BER_BP,'k-*',SNR,PER_BP,'k-+');