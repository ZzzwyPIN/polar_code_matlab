function data_crc = crccheck_m(data,gen)
[M,~] = size(data);
gen = gen(find(gen,1):end);
glen = length(gen);
dlen = length(data);
data_crc = zeros(M,glen-1);
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
    data_crc(i,:) = cr;
end
end