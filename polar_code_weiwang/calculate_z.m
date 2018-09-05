% given a z (the first entry is given), recursively calculate the other 
% entris of z

function y = calculate_z(z)
    for i = 1:log2(length(z))
        z_pre = z;
        for j=1:2^i/2
            z(2*j-1) = 2*z_pre(j)-z_pre(j)^2;
            z(2*j) = z_pre(j)^2;
        end
    end
    y=z;
end