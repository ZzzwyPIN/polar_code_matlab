%%%%%%%%%%%%%CRC ADD FUNCTION%%%%%%%%%%%%%%%
%%%This CRC add function is the plus version of the crcadd function. It
%%%can compute not only a CRC attachment of a vector but also that of a
%%%matrx.
function data_crc = crcadd_m(data,gen)
[M, N] = size(data);
gen = gen(find(gen,1):end);
glen = length(gen);
data = [data zeros(M, glen-1)];
dlen = N+glen-1;
data_crc = zeros(M,dlen);
for i = 1:M   
    cr = data(i,1:glen-1); 
    for p = glen:dlen  
        cr(glen) = data(i,p);  
        if cr(1)          
            cr = xor(cr(2:glen), gen(2:glen));   %异或运算，相当于除法     
        else
            cr = cr(2:glen);     
        end
    end
    data_crc(i,:) = [data(i,1:N) cr];
end

end