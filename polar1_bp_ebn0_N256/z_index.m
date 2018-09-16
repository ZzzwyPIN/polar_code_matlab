function index_set = z_index(n_level)
index_set = zeros(2^n_level, n_level);
    %each column
    for j = 1:n_level
       %space between each index pair
         n_sp = 2^(n_level - j);
       % groups of sections
        n_group = 2^(j-1);
       %each row
        for jj=1:n_group 
            for i = 1:n_sp
                index_set((jj-1) * n_sp * 2 + 2 * (i-1)+1,j) = (jj-1) * n_sp * 2 + i;
                index_set((jj-1) * n_sp * 2 + 2 * (i-1)+2,j) = (jj-1) * n_sp * 2 + i + n_sp;
            end
        end
    end
end

