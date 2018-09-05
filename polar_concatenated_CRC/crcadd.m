function cr=crcadd(data,gen)
glen = length(gen); 
while gen(1) == 0      
    gen = gen(2:glen);     
    glen = length(gen); 
end
data_back = data; 
for m = 1:glen-1     
    data = [data 0]; 
end
dlen = length(data);   
cr = data(1:glen-1); 
for p = glen:dlen  
    cr(glen) = data(p);  
    if cr(1)          
        cr = xor(cr(2:glen), gen(2:glen));   %异或运算，相当于除法     
    else
        cr = cr(2:glen);     
    end
end
cr = [data_back cr];
end