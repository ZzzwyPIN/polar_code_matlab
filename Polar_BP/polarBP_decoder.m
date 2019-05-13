function receive_bp_bits = polarBP_decoder(n,lr_u,lr_x,max_iter,info_index)
%%% Polar BP Decoder Function
%%% lr_x:init LLR received
%%% lr_u:init LLR send
%%% max_iter:iteration number
N = 2^n;
rel_mat_LtoR = zeros(N,n); % save polar1 decoding dataset
for idx = 1:max_iter
    % LR from Right to Left like a flood first
    rel_mat_RtoL = polar_bp_RtoL(lr_u, lr_x, rel_mat_LtoR, n);
    % Now LR from Left to right
    % Note that intrinsic information can't send to right again!
    rel_mat_LtoR = polar_bp_LtoR(lr_x, lr_u, rel_mat_RtoL, n);
end
temp = rel_mat_RtoL(:,n) + lr_u';
final_lr = temp(reverse_index(n,1:length(temp)));
receive_bp_bits = (final_lr(info_index) < 0)';