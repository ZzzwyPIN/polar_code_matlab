function z_index = get_z_index(n)
%%% this function get the index of the check nodes in each col
N = 2^n;
z_index = zeros(N,n);
for j = 1:n
    %space between each index pair
    n_sp = 2^(n - j);
    % groups of sections
    n_group = 2^(j-1);
    %each row
    for jj=1:n_group 
        for i = 1:n_sp
            z_index((jj-1) * n_sp * 2 + 2 * (i-1)+1,j) = (jj-1) * n_sp * 2 + i;
            z_index((jj-1) * n_sp * 2 + 2 * (i-1)+2,j) = (jj-1) * n_sp * 2 + i + n_sp;
        end
    end
end
    
end