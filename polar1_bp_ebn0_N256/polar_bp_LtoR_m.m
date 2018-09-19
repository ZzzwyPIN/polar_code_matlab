function rel_mat_LtoR = polar_bp_LtoR_m(rel_mat_RtoL,lr_u,lr_x,n)
N = 2^n;
rel_mat_LtoR = zeros(N,n);
z_idx = z_index(n);
% the left to right
for j = 1:n
    % the frist level
    if j == 1
         for i = 1:N/2
            % Update criterion 1
            %diagonal update top -> lower
            lr = [lr_u(z_idx(2*(i-1)+1,j)), rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1)];
            diag = checkNodeProb_sum(lr,1);
            %update lower edge
            rel_mat_LtoR(z_idx(2*(i-1)+2,j),j) = diag + lr_u(z_idx(2*(i-1)+2,j));
            % %diag node update: lower -> top: diagonal update 2
            diag = lr_u(z_idx(2*(i-1)+2,j))  + rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1);
            lr = [lr_u(z_idx(2*(i-1)+1,j)), diag];
            tmp = checkNodeProb_sum(lr,1);
            %top edge update
            rel_mat_LtoR(z_idx(2*(i-1)+1,j),j) = tmp;
         end
    else
        % the last level
        if j==n
            for i=1:N/2
%                     diagonal update top -> lower
                lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), lr_x(z_idx(2*(i-1)+1,j))];
                diag = checkNodeProb_sum(lr,1);
%                     update lower edge
                rel_mat_LtoR(z_idx(2*(i-1)+2,j),j) = diag + rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1);
%                     diag node update: lower -> top: diagonal update 2
                diag =  rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1) + lr_x(z_idx(2*(i-1)+2,j));
                lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), diag];
                tmp = checkNodeProb_sum(lr,1);
%                     update top edge
                rel_mat_LtoR(z_idx(2*(i-1)+1,j),j) = tmp;
            end
        else
            % intermediate level
            for i=1:N/2
                %diagonal update top -> lower
                lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), rel_mat_RtoL(z_idx(2*(i-1)+1,j),j+1)];
                diag = checkNodeProb_sum(lr,1);
                %store lower edge 
%                             bot_old = rel_mat_LtoR(z_idx(2*(i-1)+2,j),j);
                %update lower edge
                rel_mat_LtoR(z_idx(2*(i-1)+2,j),j) = diag + rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1);
                % diag node update: lower -> top: diagonal update 2
                diag =  rel_mat_LtoR(z_idx(2*(i-1)+2,j),j-1) + rel_mat_RtoL(z_idx(2*(i-1)+2,j),j+1);
                lr = [rel_mat_LtoR(z_idx(2*(i-1)+1,j),j-1), diag];
                tmp = checkNodeProb_sum(lr,1);
                rel_mat_LtoR(z_idx(2*(i-1)+1,j),j) = tmp;
            end
        end
    end
end