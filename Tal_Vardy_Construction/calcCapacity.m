function c = calcCapacity( a, b )
% calculate capacity given by a and b
%   c = -(a+b)log2((a+b)/2) + a log2(a) + b log2(b)
  
if (a == 0) 
    if b == 0
        c = 0;
    else
        c = -(b)*log2((b)/2) + b*log2(b);
    end
else
    if b == 0
        c = -(a)*log2((a)/2) + a*log2(a);
    else
        c = -(a+b)*log2((a+b)/2) + a*log2(a) + b*log2(b);
    end
end

end

