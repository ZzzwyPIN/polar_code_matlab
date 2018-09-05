function s_temp2 = insert_bit(s_temp1,s_temp2,temp_index1,temp_index2)
%%%%%Insert the bad information bits in Polar1 to the good in polar2.
%   s_temp1: bitsource of polar1
%   s_temp2: bitsource of polar2 before compeleting
%   temp_index1: The index of s_temp1 to insert to s_temp2
%   temp_index2: The index of s_temp2 to recive inserted bits
for i = 1:length(temp_index1)
    s_temp2(temp_index2(i)+1:length(s_temp2)+1) = s_temp2(temp_index2(i):end);
    s_temp2(temp_index2(i)) = s_temp1(temp_index1(i));
end