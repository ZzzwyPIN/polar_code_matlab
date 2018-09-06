clc
clear

% ������������
n = 8;  % ����λ��
R = 0.2188;    % ����

SNR = -3:2;

init_lr_max = 3;    % limit the max LR of the channel to be with [-3 3]
max_iter = 30;
block_num = 10000;

% ��������
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
Info_index = sort(I(K:-1:1));  % ��ѡ�����õ��ŵ�������Ϣλ
Frozen_index = sort(I(end:-1:K+1));   % ���䶳��λ���ŵ�

% get generate matrix
G = encoding_matrix(n);
Gi = G(Info_index,:);
Gf = G(Frozen_index,:);
frozen_bits = zeros(1,k_f);
rng('shuffle')
for i = 1:length(SNR)
    sigma = (1/snr(i))^0.5;
    % set PER and BER counter
    PerNum = 0;
    BerNum = 0;
    for iter = 1:block_num
        fprintf('\nNow iter: %2d\tNow SNR: %d', iter, SNR(i));
        source_bit = randi([0 1],1,K);
        encode_temp = rem(source_bit*Gi + frozen_bits*Gf,2);
    
        % bpsk modulation
        encode_temp = (-1).^(encode_temp + 1);
        % add noise
        receive_sample = encode_temp + sigma * randn(size(encode_temp));
        
        lr_x = -2*receive_sample./(sigma^2);
        % decoding follow
        lr_u = zeros(1,N); % save send sample LR in each iteration
        lr_u(reverse_index(n,Frozen_index)) = init_max;
        
        receive_bits = polarBP_decoder(n,lr_u,lr_x,max_iter,Info_index);
        
        % calculate BER and PER
        count = sum(receive_bits ~= source_bit);
        if count ~= 0
            PerNum = PerNum + 1;
            BerNum = BerNum + count;
        end 
    end
    perBP(i) = PerNum/block_num;
    berBP(i) = BerNum/(K*block_num);
end