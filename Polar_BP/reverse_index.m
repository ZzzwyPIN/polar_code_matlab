function index = reverse_index(n,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is used to compute the reversation of bit source's index
%%% and it is the index of the channel input index after encoding.
[M, N] = size(x);
if (M > 1 && N > 1)
  error(message('x is an Array'));  
end
if (M == 1 && N > 1)
    index = zeros(size(x));
    for i = 1:N
        temp = x(i) - 1;
        idx = dec2bin(temp);
        n_idx = length(idx);
        for k=1:n-n_idx
            idx = strcat('0',idx);
        end
        index(i) = bin2dec(fliplr(idx))+1;
    end
end
if (M >1 && N == 1)
    index = zeros(size(x));
    for i = 1:M
        temp = x(i) - 1;
        idx = dec2bin(temp);
        n_idx = length(idx);
        for k=1:n-n_idx
            idx = strcat('0',idx);
        end
        index(i) = bin2dec(fliplr(idx))+1;
    end
end
if (M == 1 && N == 1)
    temp = x - 1;
    idx = dec2bin(temp);
    n_idx = length(idx);
    for k=1:n-n_idx
        idx = strcat('0',idx);
    end
    index = bin2dec(fliplr(idx))+1;
end
