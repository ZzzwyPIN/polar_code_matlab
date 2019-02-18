%%% LL June 11, 2015
% calculate the error probability for each bit channel
% sort the bit channels according to their error probability
% select the best bit channels as the information set

%% Aaron Wong July, 2015
% dealing with continuous AWGN channel

%%
clc;clear;
n=8; % coding level
N = 2^n; % block length

CHANNEL_TYPE = 0; % 0: AWGN, 1: BSC, 2: BEC

R = 0.75; % code rate 
K = floor(N * R); % # of inforamtion bits in one code block
code_idx = zeros(1,K); % # of frozen bits in one code block
frozen_idx = zeros(1,N-K);

TRAN_COMBINE = 1;
log_en = 0;

% # of output alphabet set, 
v = 2^6;
u = 2*v; % final channel output alphabet size
% size of the degrading of the continuous AWGN channel to a discrete channel with a
% finite output size u=2*v 
v_awgn = 2^10;
u_awgn = 2 * v_awgn;

% BSC channel
if CHANNEL_TYPE == 1
    p = 0.1; % transition probability
    % underlying channe W: BSC
    L = 1;
    W = [1-p, p]; % W stores W(y|0) W(\bar(y)|0)
    % we don't need the LR info as the input
    %LR = (1-p)/p; %  stores the sorted LR
end
%%
% this is degrade specific!
% AWGN channel
% use u as output size so that difference to original channel is bounded by 2/u
%yl: lower bound of y in a symbol
%yh: higher bound of y in a symbol
if CHANNEL_TYPE == 0  
   %param
   SNR = 4.6494; % in dB
   snr = 10^(SNR/10);
   esn0 = snr * R;
   sigma = (2*esn0)^-0.5;
%    sigma = snr^-0.5;
%    sigma = 0.1;
   % size of the degrading of the continuous AWGN channel to a discrete channel with a
   % finite output size u=2*v
%    v = 2^6;
%    u=2*v;
   W = zeros(1,u_awgn);
   %LR = zeros(1,v);
   
   %since we know 0 == calcC(0,sigma) is true
   yl = 0;
   yh = zeros(1,v_awgn);
   yh(1) = yl;
   %loop
   lamda_next = 1;
   for i=1:v_awgn-1
        %now we need to solve
        % i/v == calcC(y,sigma)
        % to determine of y of this symbol.
        % but the function is nonlinear!
        % suppose such function exists:(yl is for optimization)
        % lamda_next is returned to be used as the starting point for the
        % next round
        [yh(i+1),lamda_next] = revcalcC(i/v_awgn,lamda_next,sigma);
        %then we get
        % the area under a range is obtained from the CDF, not PDF
%         wy0 = normpdf(yh,0,sigma) - normpdf(yl,0,sigma);
%         wy1 = normpdf(yh,1,sigma) - normpdf(yl,1,sigma);
        wy0 = normcdf(yh(i+1),1,sigma) - normcdf(yl,1,sigma);
        wy1 = normcdf(yh(i+1),-1,sigma) - normcdf(yl,-1,sigma);
        %we place it in such an order
 %       W[2*(i-1)] = wy1;
%        w[2*(i-1)+1] = wy0;
        % to facilitate the merging, the transition probabilities are
        % placed in this way, y1,y2,...,yL, y1c,y2c,...,yLc
        W(i) = wy0;
        W(i+v_awgn) = wy1;
        %LR(i) = wy0/wy1;
        %prepare for next symbol
        yl = yh(i+1);
   end
   
   %last symbol calculation
   %we y will go to infinite at the last symbol so
%    wy0 = 1 - normpdf(yl,0,sigma);
%    wy1 = 1 - normpdf(yl,1,sigma);
   wy0 = 1 - normcdf(yl,1,sigma);
   wy1 = 1 - normcdf(yl,-1,sigma);
   %we place it in such an order
%    W[2*(v-1)] = wy1;
%    w[2*(v-1)+1] = wy0;
   W(v_awgn) = wy0;
   W(2*v_awgn) = wy1;
   %LR(v) = wy0/wy1;
end
%%
%filename = ['awgn_snr0p1875db_',num2str(u_awgn)];
%save(filename,'W');

% given a bit channel i
%for i=1:N
for i=1:N
    fprintf('\nNow iter : %3d',i);
    % Initial processing 
    W1 = sortTran_sim(W);
    [Q] = degrading_merge(W1, u);
    bi = dec2bin(i-1,n);
    for m=1:n 
        nw = length(Q);
        % sort the LRs according to bi(m)
        if bi(m) == '0'
            % size of the output alphabet at this level 
            Nm = nw^2;
            if TRAN_COMBINE
                [W10] = calcTran_sim_combine(Q,0); % type 0 transformation 
            else
                [W10] = calcTran_sim(Q,0); % type 0 transformation 
            end
            % calculate LRs for W
            nw = length(W10);
            L = nw/2;
            [W20] = sortTran_sim(W10);
            [Q] = degrading_merge(W20,u);
        else
            Nm = nw^2*2;
            if TRAN_COMBINE
                [W11] = calcTran_sim_combine(Q,1); % type 1 transformation 
            else
                [W11] = calcTran_sim(Q,1); % type 1 transformation 
            end
            % calculate LRs for W
            nw = length(W11);
            L = nw/2;
            [W21] = sortTran_sim(W11);
            [Q] = degrading_merge(W21,u);
        end
    end
    nq = length(Q);
    P(i) = sum(Q(nq/2+1:nq))/2;
    
    
    
    % calculate the mutual infomation
    % w(y)
    Wy = zeros(1,nq/2);
    Wy = (Q(1:nq/2) + Q(nq/2+1:nq))/2;
    % Hy
    Hy = -2*sum(Wy .* log2(Wy));
    % H(Y|X)
    Hyx = -sum(Q.*log2(Q));
    Ixy(i) = Hy-Hyx;
    
    if log_en
        path = './log_files/';
        filename = [path, 'bit_channel_n_', num2str(N),'_i_',num2str(i)];
        if TRAN_COMBINE
            filename = [filename,'_comb'];
        end
        filename = [filename, '.txt'];
        dlmwrite(filename,Q,' ');
    end
end
%%
% [P_s,I] = sort(P);
% code_idx = I(1:K);
% frozen_idx = I(K+1:N);
% code_idx = sort(code_idx);
% frozen_idx = sort(frozen_idx);
% 
% if log_en
%     filename = [path,'bit_channel_Pe_', num2str(N),];
%     if TRAN_COMBINE
%         filename = [filename,'_comb'];
%     end
%     filename = [filename, '.txt'];
%     dlmwrite(filename,P,' ');
% end
% 
% path = './log_files/';
% filename = [path,'code_idx_D_AWGNC0dB_Rr0p5_n_', num2str(n),'.dat'];
% dlmwrite(filename,code_idx')
% filename = [path,'frozen_idx_D_AWGNC0dB_Rr0p5_n_', num2str(n),'.dat'];
% dlmwrite(filename,frozen_idx')





