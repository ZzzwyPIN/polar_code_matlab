function rel_mat_RtoL = polar_bp_RtoL(lr_u, lr_x, rel_mat_LtoR, n)
%%% this function decode the polar code with BP
%   Frozen_index: the index of frozen bits
%   lr_u: the LR of send sample
%   n: bit order
N = 2^n;
rel_mat_RtoL = zeros(N,n);
z_index = get_z_index(n);
% each col
for j = 1:n
    if j == n
        for i = 1:N/2
            % low to top
            diag = lr_u(z_index(2*(i-1)+2,1)) + lr_x(z_index(2*(i-1)+2,1));
            lr = [lr_x(z_index(2*(i-1)+1,1)) diag];
            rel_mat_RtoL(z_index(2*(i-1)+1,1),j) = checkNodeProbSum(lr,1);
            % top to low
            lr = [lr_x(z_index(2*(i-1)+1,1)) lr_u(z_index(2*(i-1)+1,1))];
            diag = checkNodeProbSum(lr,1);
            rel_mat_RtoL(z_index(2*(i-1)+2,1),j) = lr_x(z_index(2*(i-1)+2,1)) + diag;
        end
    else
        for i = 1:N/2
            % low to top   
            diag = rel_mat_LtoR(z_index(2*(i-1)+2,n-j+1),n-j) + lr_x(z_index(2*(i-1)+2,n-j+1));
            lr = [lr_x(z_index(2*(i-1)+1,n-j+1)) diag];
            rel_mat_RtoL(z_index(2*(i-1)+1,n-j+1),j) = checkNodeProbSum(lr,1);
            % top to low
            lr = [lr_x(z_index(2*(i-1)+1,n-j+1)) rel_mat_LtoR(z_index(2*(i-1)+1,n-j+1),n-j)];
            diag = checkNodeProbSum(lr,1);
            rel_mat_RtoL(z_index(2*(i-1)+2,n-j+1),j) = lr_x(z_index(2*(i-1)+2,n-j+1)) + diag;
        end
        lr_x = rel_mat_RtoL(:,j);
    end
end