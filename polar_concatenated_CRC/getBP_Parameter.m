function [lr_u,lr_x] = getBP_Parameter(receive_sample,frozen_bits,frozen_index,n,init_max,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%this founction receive the channel outputs and frozen_index
%%%and return the LLR and init lr_u for BP decoding
N = 2^n;
lr_u = zeros(1,N); % save send sample LR in each iteration
frozen_index_0 = frozen_bits == 0;
frozen_index_1 = frozen_bits == 1;
lr_u(reverse_index(n,frozen_index(frozen_index_0))) = init_max;
lr_u(reverse_index(n,frozen_index(frozen_index_1))) = -init_max;
lr_x = -2*receive_sample./(sigma^2);