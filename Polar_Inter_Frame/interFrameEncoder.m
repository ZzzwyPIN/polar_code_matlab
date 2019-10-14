function encoder_temp = interFrameEncoder(u, lambda_offset, llr_layer_vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,M] = size(u);
encoder_temp = zeros(N,M);
for i = 1:M
    encoder_temp(:,i) = polar_encoder(u(:,i), lambda_offset, llr_layer_vec);
end
