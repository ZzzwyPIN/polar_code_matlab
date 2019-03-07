function s_temp2 = insert_bit(s_temp1,s_temp2,temp_index1,temp_index2)
%%%%%nargin = 3: insert vector s_temp1 into s_temp2 in an order of
%%%%%temp_index1
%%%%%nargin = 4:Insert the bad information bits in Polar1 to the good in polar2.
%   s_temp1: bitsource of polar1
%   s_temp2: bitsource of polar2 before compeleting
%   temp_index1: The index of s_temp1 to insert to s_temp2
%   temp_index2: The index of s_temp2 to recive inserted bits
if nargin == 3
    for i = 1:length(temp_index1)
        s_temp2(temp_index1(i)+1:length(s_temp2)+1) = s_temp2(temp_index1(i):end);
        s_temp2(temp_index1(i)) = s_temp1(i);
    end
else
    for i = 1:length(temp_index1)
        s_temp2(temp_index2(i)+1:length(s_temp2)+1) = s_temp2(temp_index2(i):end);
        s_temp2(temp_index2(i)) = s_temp1(temp_index1(i));
    end
end