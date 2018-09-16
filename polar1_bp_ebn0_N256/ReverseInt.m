function out = ReverseInt(x,x_len )
    
    out = 0;
    for i=0:(x_len-1)
        %out = out | (((x & (1 << i)) >> i) << (x_len - 1 - i));
        out = bitor(out,bitshift((bitshift((bitand(bitshift(1,i),x)),-i)),(x_len - 1 - i)));
    end
    
end

