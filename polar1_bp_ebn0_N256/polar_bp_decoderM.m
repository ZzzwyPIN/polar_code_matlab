function decision_bits = polar_bp_decoderM(received_sample,sigma,info_index,frozen_bits,frozen_index,max_iter,init_max,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BP decoder
N = 2^n;
init_lr = -2 * received_sample/(sigma^2);
 % establish the check nodes / variable nodes connections
rel_mat_RtoL = zeros(N,n);
rel_mat_LtoR = zeros(N,n);
lr_u = zeros(1,N); % save send sample LR in each iteration
frozen_index_0 = frozen_bits == 0;
frozen_index_1 = frozen_bits == 1;
lr_u(reverse_index(n,frozen_index(frozen_index_0))) = init_max;
lr_u(reverse_index(n,frozen_index(frozen_index_1))) = -init_max;

% LR values fed to the graph from the right: initial LRs from the channel
lr_x = init_lr;
% LR values fed to the graph from the right, used when calculating
for n_i=1:max_iter
    rel_mat_RtoL = polar_bp_RtoL_m(rel_mat_LtoR,lr_x,lr_u,n);
    rel_mat_LtoR = polar_bp_LtoR_m(rel_mat_RtoL,lr_u,lr_x,n);
end
temp = rel_mat_RtoL(:,1)+init_lr';
% received bits
final_lr = temp(reverse_index(n,1:length(temp)));
decision_bits = (final_lr(info_index) < 0)';
