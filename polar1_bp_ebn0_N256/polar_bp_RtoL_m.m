function rel_mat_RtoL = polar_bp_RtoL_m(rel_mat_LtoR,lr_x,lr_u,n)
N = 2^n;
rel_mat_RtoL = zeros(N,n);
z_idx = z_index(n);
% Right to left    
for j = n:-1:1
% the first level
    if j==1
        for i=1:N/2
            % top node update
            % diagonal update: lower -> top
            diag = rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1) + lr_u(z_idx(2*(i-1)+2,j));
            lr = [rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1), diag];
            % top node update
            rel_mat_RtoL(z_idx(2*(i-1)+1,j),j) = checkNodeProb_sum(lr,1);
            % bottom node update
            % diagonal update: top -> lower
            lr = [lr_u(z_idx(2*(i-1)+1,j)), rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1)];
            diag = checkNodeProb_sum(lr,1); 
             % bottom edge update
            rel_mat_RtoL(z_idx(2*(i-1)+2,j),j) = rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1)+diag;
        end
    else
        % the last level
       if j==n
            for i=1:N/2
                % diagonal update: lower -> top
                diag = lr_x(z_idx(2*(i-1)+2,j))+rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1);
                lr = [lr_x(z_idx(2*(i-1)+1,j)), diag];
                 % top edge update
                rel_mat_RtoL(z_idx(2*(i-1)+1,j),j) = checkNodeProb_sum(lr,1);
                % diagonal update: top -> lower
                lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), lr_x(z_idx(2*(i-1)+1,j))];
                diag = checkNodeProb_sum(lr,1); 
                 % bottom edge update
                rel_mat_RtoL(z_idx(2*(i-1)+2,j),j) = lr_x(z_idx(2*(i-1)+2,j))+diag;
            end
       else
           % intermediate level
           for i=1:N/2
                % diagonal update: lower -> top
                diag = rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1)+rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1);
                lr = [rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1), diag];
                 % top edge update
                rel_mat_RtoL(z_idx(2*(i-1)+1,j),j) = checkNodeProb_sum(lr,1);
                % diagonal update: top -> lower
                lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1)];
                diag = checkNodeProb_sum(lr,1); 
                 % bottom edge update
                rel_mat_RtoL(z_idx(2*(i-1)+2,j),j) = rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1)+diag;
            end
       end
    end
end