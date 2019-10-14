function u = getReverseBit(N, M, pure_info_index, MUUB, MRFB, poly, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(N,M);

info_index = [MRFB MUUB pure_info_index];

infoNum = length(info_index);
r = length(poly)-1;
source_bits = rand(infoNum - m - r,1)>0.5;
source_bits_crc = crcadd(source_bits,poly);
u(info_index, 1) = source_bits_crc;

for i = 2:M
    source_bits = rand(infoNum - m - r,1)>0.5;
    source_bits_crc = crcadd(source_bits,poly);
    u(info_index, i) = source_bits_crc;
    u(MRFB,i) = u(MUUB,i-1);
end


