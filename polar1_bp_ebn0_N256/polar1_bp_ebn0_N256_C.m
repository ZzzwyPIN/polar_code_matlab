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
g_bp_max = g_bp_init_max * (n); % used to indicate max LR of frozen bits
%g_bp_inter_max = g_bp_max + n*g_bp_init_max; % intermediate max LR
if (g_bp_max > 30)
   g_bp_max = 30;
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
        init_lr = -2 * received_sample/(sigma^2);

        % establish the check nodes / variable nodes connections
        z_idx = z_index(n);


        % used to store the LRs from right to left and from left to right
        rel_mat_RtoL = zeros(N,n);
        rel_mat_LtoR = zeros(N,n);
        init_frozen_lr = zeros(1,N);
        for i=1:length(frozen_bits)
            if(frozen_bits(i) == 0)
                init_frozen_lr(ReverseInt(frozen_index(i)-1,n)+1) = g_bp_max;
            else
                init_frozen_lr(ReverseInt(frozen_index(i)-1,n)+1) = -g_bp_max;
            end
        end

        % calculates all bit indices in bit reversed order
        code_idx_rev = zeros(1,N);
        for i=1:N
            code_idx_rev(i) = ReverseInt(i-1,n) +1;
        end
        % LR values fed to the graph from the left: only contains LR values from
        % frozen bits
        lr_u = init_frozen_lr';
        % LR values fed to the graph from the right: initial LRs from the channel
        lr_x = init_lr';
        % LR values fed to the graph from the right, used when calculating
        % rel_mat_LtoR. In the simulation, it's made the same as lr_x
        lr_x_LtoR = lr_x;
        % LR values fed to the graph from the left, used when calculating
        % rel_mat_LtoR. In the simulation, it's made the same as lr_u
        lr_u_LtoR = lr_u;
        for n_i=1:max_iter
            % Right to left    
            for j=n:-1:1
            % the first level
                if j==1
                    for i=1:N/2
                        % top node update
                        % diagonal update: lower -> top
                        diag = rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1) + lr_u(z_idx(2*(i-1)+2,j));
                        lr = [rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1), diag];
                        % top node update
                        rel_mat_RtoL(z_idx(2*(i-1)+1,j),j) = checkNodeProb_sum(lr,1);
                        % bottom node update
                        % diagonal update: top -> lower
                        lr = [lr_u(z_idx(2*(i-1)+1,j)), rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1)];
                        diag = checkNodeProb_sum(lr,1); 
                         % bottom edge update
                        rel_mat_RtoL(z_idx(2*(i-1)+2,j),j) = rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1)+diag;
                    end
                else
                    % the last level
                   if j==n
                        for i=1:N/2
                            % diagonal update: lower -> top
                            diag = lr_x(z_idx(2*(i-1)+2,j))+rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1);
                            lr = [lr_x(z_idx(2*(i-1)+1,j)), diag];
                             % top edge update
                            rel_mat_RtoL(z_idx(2*(i-1)+1,j),j) = checkNodeProb_sum(lr,1);
                            % diagonal update: top -> lower
                            lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), lr_x(z_idx(2*(i-1)+1,j))];
                            diag = checkNodeProb_sum(lr,1); 
                             % bottom edge update
                            rel_mat_RtoL(z_idx(2*(i-1)+2,j),j) = lr_x(z_idx(2*(i-1)+2,j))+diag;
                        end
                   else
                       % intermediate level
                       for i=1:N/2
                            nodes_id = [z_idx(2*(i-1)+1,j),z_idx(2*(i-1)+2,j)];
                            % diagonal update: lower -> top
                            diag = rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1)+rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1);
                            lr = [rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1), diag];
                             % top edge update
                            rel_mat_RtoL(z_idx(2*(i-1)+1,j),j) = checkNodeProb_sum(lr,1);
                            % diagonal update: top -> lower
                            lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1)];
                            diag = checkNodeProb_sum(lr,1); 
                             % bottom edge update
                            rel_mat_RtoL(z_idx(2*(i-1)+2,j),j) = rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1)+diag;
                        end
                   end
                end
            end
            
            
            % the left to right
            for j=1:n
                % the frist level
                if j==1
                     for i=1:N/2
                        % Update criterion 1
                        %diagonal update top -> lower
                        lr = [lr_u_LtoR(z_idx(2*(i-1)+1,j)), rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1)];
                        diag = checkNodeProb_sum(lr,1);
                        %update lower edge
                        rel_mat_LtoR(z_idx(2*(i-1)+2,j),j) = diag + lr_u_LtoR(z_idx(2*(i-1)+2,j));
                        % %diag node update: lower -> top: diagonal update 2
                        diag = lr_u_LtoR(z_idx(2*(i-1)+2,j))  + rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1);
                        lr = [lr_u_LtoR(z_idx(2*(i-1)+1,j)), diag];
                        tmp = checkNodeProb_sum(lr,1);
                        %top edge update
                        rel_mat_LtoR(z_idx(2*(i-1)+1,j),j) = tmp;
                     end
                else
                    % the last level
                    if j==n
                        for i=1:N/2
        %                     diagonal update top -> lower
                            lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), lr_x_LtoR(z_idx(2*(i-1)+1,j))];
                            diag = checkNodeProb_sum(lr,1);
        %                     update lower edge
                            rel_mat_LtoR(z_idx(2*(i-1)+2,j),j) = diag + rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1);
        %                     diag node update: lower -> top: diagonal update 2
                            diag =  rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1) + lr_x_LtoR(z_idx(2*(i-1)+2,j));
                            lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), diag];
                            tmp = checkNodeProb_sum(lr,1);
        %                     update top edge
                            rel_mat_LtoR(z_idx(2*(i-1)+1,j),j) = tmp;
                        end
                    else
                        % intermediate level
                        for i=1:N/2
                            %diagonal update top -> lower
                            lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1)];
                            diag = checkNodeProb_sum(lr,1);
                            %store lower edge 
%                             bot_old = rel_mat_LtoR(z_idx(2*(i-1)+2,j),j);
                            %update lower edge
                            rel_mat_LtoR(z_idx(2*(i-1)+2,j),j) = diag + rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1);
                            % diag node update: lower -> top: diagonal update 2
                            diag =  rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1) + rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1);
                            lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), diag];
                            tmp = checkNodeProb_sum(lr,1);
                            rel_mat_LtoR(z_idx(2*(i-1)+1,j),j) = tmp;
                        end
                    end
                end
            end
        end
        % LRs used for decisions
        temp = rel_mat_RtoL(:,1)+init_frozen_lr';
        % received bits
        final_lr = temp(reverse_index(n,1:length(temp)));
        receive_bits = (final_lr(info_index) < 0)';
        
        count = sum(receive_bits ~= source_bits);
        if count ~= 0
           ber = ber + count;
           per = per +1;
        end
    end
    BER_BP(1,i_aw)= ber/(K * block_num);
    PER_BP(1,i_aw)= per / block_num;
    disp('循环次数')
    i_aw
end
BER_BP
PER_BP


