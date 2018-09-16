% BP decoding of polar codes
% Dr Jan. 6, 2017
clc; clear;

% diary on
% TOP_SCALE = 0;
n = 8;
N = 2^n;
R = 0.36; 
block_num = 10000;
SNR=0;
snr=10.^(SNR/10);
g_bp_init_max = 3; % limit the max LR of the channel to be within [-3, 3]
g_bp_max = g_bp_init_max * (n); % used to indicate max LR of frozen bits
%g_bp_inter_max = g_bp_max + n*g_bp_init_max; % intermediate max LR
if (g_bp_max > 30)
   g_bp_max = 30;
end
max_iter=30;
k=floor(N*R);%infromation bits in one block
k_f=N-k;%frozen_bits in one block
G= encoding_matrix(n);              %产生G矩阵   F@(i)F     (1;n+1;2;n+2;3;n+3;......;n-1;2n-1;n;2n)


% source_bits=zeros(1,k);
frozen_bits=randi([0,1],1,N-k);
encoded_bits=zeros(1,N);%极化码编码矩阵
received_bits=zeros(1,N);
% load('Pe_snr1p0db_R0p25_n_322048n_5.mat');
load('Pe_snr0p0db_2048_n_8.mat');
Ptmp = sort(P); 
% (2) find cut-off Z value below which the rows of G will be used for encoding
cut_p = Ptmp(k); %Ztmp(4)
% (3) find row number of G used for encoding 
coding_idx = find(P <= cut_p);   %编码 Z<=cut_z  传输速率小的信道容量大，信道好
% (4) find the frozen bit index  
frozen_idx = find(P > cut_p);
Gi = G(coding_idx,:);  % 选取G矩阵coding行数
% (6) form the encoding matrix for frozen bits
Gf = G(frozen_idx,:);  % 选取G矩阵frozen行数
rng('shuffle');
ebn0=k/N*snr;
for i_aw=1:length(SNR)
    sigma=(2*ebn0(i_aw))^(-0.5);
    ber = 0;
    per = 0;
%     for j=1:length(source_bits)/k  %1：（40/4）
%         encoded_bits(1,(j-1)*N+1:j*N) = source_bits((j-1)*k+1:j*k)*Gi...   
%             + frozen_bits((j-1)*k_f+1:j*k_f)*Gf;    %encoded_bits=source_bits*Gi+frozen_bits*Gf   
%     end
    % BPSK modulation
    
    for r_idx = 1:block_num
        fprintf('\nNow iter: %d\tNow SNR: %d',r_idx, SNR(i_aw));
%         decoded_info_one=zeros(1,k);
        init_lr=zeros(1,N);
        noise=sigma*randn(1,N);
        source_bits=randi([0,1],1,k);
        encoded_bits = source_bits*Gi+frozen_bits*Gf;
        idx_0 = find(mod(encoded_bits,2) ==  0);
        idx_1 = find(mod(encoded_bits,2) ==  1);
        encoded_bits(idx_0) = -1*ones(1,length(idx_0));
        encoded_bits(idx_1) = ones(1,length(idx_1));
        received_sample = encoded_bits + noise;%noise is added
        for j=1:N  
            init_lr(j)=exp(-2*received_sample(j)*(2*ebn0(i_aw)));
            tmp = init_lr(j);
            if (tmp == 0)
%                 error('value error!');
                init_lr(j) = -1 * g_bp_max;
            else 
                init_lr(j) = log(tmp);
            end 
            if (init_lr(j) < -30) 
                init_lr(j) = -1 * g_bp_max;
            end
        end    
        %the following limit will cause the bad performance of BP decoding in large SNR area.
%          Initial LR preprocessing to make it within [-3, 3]
%         max_init = max(abs(init_lr));
%         if max_init > g_bp_init_max
%             init_lr = init_lr / max_init * g_bp_init_max;
%         end

        % establish the check nodes / variable nodes connections
        z_idx = z_index(n);

        % calculate the code index and frozen index according to the Z values
        % [coding_idx,frozen_idx] = coding_index(0.1,n,R)

        % used to store the LRs from right to left and from left to right
        rel_mat_RtoL = zeros(N,n);
        rel_mat_LtoR = zeros(N,n);

        % used to store the received bits in nature order
        received_bits_one = zeros(1,N);
        % used to store the frozen indices in the bit reversed order
        frozen_idx_rev = zeros(1,length(frozen_idx));

        % used to store the first level of LRs known frozen bits, 
        % the info bits LR is 0
        init_frozen_lr = zeros(1,N);
        % used to store the final LR values for decisions 
        final_lr = zeros(1,N)';
        frozen_bits_one=frozen_bits;
        % the LR info fed to the variable nodes of level 1 
        for i=1:length(frozen_bits_one)
            if(frozen_bits_one(i) == 0)
                init_frozen_lr(ReverseInt(frozen_idx(i)-1,n)+1) = g_bp_max;
            else
                init_frozen_lr(ReverseInt(frozen_idx(i)-1,n)+1) = -g_bp_max;
            end
            % feed frozen bits to received_bits
            received_bits_one(frozen_idx(i)) = frozen_bits_one(i);
            % calculates the frozen bit indices in bit reversed order
            frozen_idx_rev(i) = ReverseInt(frozen_idx(i)-1,n)+1;
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
            n_i;
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

%             if TOP_SCALE  the procedure can't reach here. it is useless
%                 rel_max = max(max(abs(rel_mat_RtoL)));
%                 if(rel_max > g_bp_max)
%                     rel_mat_RtoL = rel_mat_RtoL/ rel_max * g_bp_max;
%                     rel_mat_LtoR = rel_mat_LtoR / rel_max * g_bp_max;
%                 end
%             end
        %     disp('Right->Left')
            rel_mat_RtoL;
            % LRs used for decisions
            final_lr = rel_mat_RtoL(:,1)+init_frozen_lr';
            % received bits
            for i=1:length(coding_idx)
                idx_r = ReverseInt(coding_idx(i)-1,n) + 1;
                if(final_lr(idx_r) < 0)
                    received_bits_one(coding_idx(i)) = 1;
                else
                    received_bits_one(coding_idx(i)) = 0;
                end
            end

            % decide the code word
            x_bits = zeros(1,N);
            for i=1:N
                if((rel_mat_LtoR(i,n)+init_lr(i)) < 0)
                    x_bits(i) = 1;
                end
            end
            x_bits;

            % left to Right
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
                            bot_old = rel_mat_LtoR(z_idx(2*(i-1)+2,j),j);
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
%            if TOP_SCALE  the procedure can't reach here. it is useless
%                rel_max = max(max(abs(rel_mat_LtoR)));
%                if(rel_max > g_bp_LtoR_max)
%                    rel_mat_LtoR = rel_mat_LtoR/ rel_max * g_bp_LtoR_max;
%                    rel_mat_RtoL = rel_mat_RtoL/ rel_max * g_bp_LtoR_max;
%                end
%            end
        %    disp('Left->Right after Scale')
%            rel_mat_LtoR;
%            % lR used for input bits
%            final_lr;
           % LRs used for code word
           final_lr_x = rel_mat_LtoR(:,n) + init_lr';
        end
        for i=1:length(coding_idx)
            idx_re(i)= ReverseInt(coding_idx(i)-1,n) + 1;
        end
        final_lr_info=final_lr(idx_re);
        sinf_mem_lr(k*(r_idx-1)+1:k*(r_idx-1)+k)=final_lr_info(1:k);%get all the information bits' LR values
        decoded_info_one = received_bits_one(coding_idx);%get decoded bits of one block 
%         decoded_info((r_idx-1)*k+1:r_idx*k)=decoded_info_one;%get the decoded bits of all polar blocks
%         gp_BER(r_idx)=sum(decoded_info_one~=x(r_idx,:))/k; %get the BER of each block
%         if gp_BER(r_idx)~=0
%             err_group(r_idx)=r_idx;
%         end
        
        count = sum(decoded_info_one ~= source_bits);
        if count ~= 0
           ber = ber + count;
           per = per +1;
        end
    end
    BER_BP(1,i_aw)= ber/(k * block_num);
    PER_BP(1,i_aw)= per / block_num;
    disp('循环次数')
    i_aw
end
BER_BP
PER_BP


