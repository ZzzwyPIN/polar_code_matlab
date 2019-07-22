function rel_mat_LtoR = polar_bp_LtoR(lr_x, lr_u, rel_mat_RtoL, n)
%%% this function decode the polar code with BP 
%   rel_mat_LtoR: 
%   lr_x: the LR of receive sample
%   n: bit order
N = 2^n;
rel_mat_LtoR = zeros(N,n);
z_index = get_z_index(n);
% each col
for j = 1:n
    if j == n
        for i = 1:N/2
            % low to top
            diag = lr_u(z_index(2*(i-1)+2,j)) + lr_x(z_index(2*(i-1)+2,j));
            lr = [lr_u(z_index(2*(i-1)+1,n-j+1)) diag];
            rel_mat_LtoR(z_index(2*(i-1)+1,j),j) = checkNodeProbSum(lr,1);
            % top to low
            lr = [lr_u(z_index(2*(i-1)+1,j)) lr_x(z_index(2*(i-1)+1,j))];
            diag = checkNodeProbSum(lr,1);
            rel_mat_LtoR(z_index(2*(i-1)+2,j),j) = lr_u(z_index(2*(i-1)+2,j)) + diag;
        end
    else
        for i = 1:N/2
            % low to top   
            diag = rel_mat_RtoL(z_index(2*(i-1)+2,j),n-j) + lr_u(z_index(2*(i-1)+2,j));
            lr = [lr_u(z_index(2*(i-1)+1,j)) diag];
            rel_mat_LtoR(z_index(2*(i-1)+1,j),j) = checkNodeProbSum(lr,1);
            % top to low
            lr = [lr_u(z_index(2*(i-1)+1,j)) rel_mat_RtoL(z_index(2*(i-1)+1,j),n-j)];
            diag = checkNodeProbSum(lr,1);
            rel_mat_LtoR(z_index(2*(i-1)+2,j),j) = lr_u(z_index(2*(i-1)+2,j)) + diag;
        end
        lr_u = rel_mat_LtoR(:,j);
    end
end