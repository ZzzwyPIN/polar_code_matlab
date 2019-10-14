function [Q,conjIndex] = calcTran( W,cIndex, type_flag )
% calcTran calculate the transition probability of the one-step channel transformation
%   W is the original channel transition probabilities
% type_flag: 0, type 0 transform; 1: type 1 transform

L = length(W)/2;

if type_flag == 0 % type 0 transformation, Q's size squares that of W
    nq = (2*L)^2;
    done_flag = zeros(1,nq);
    Q = zeros(1,nq);
    conjIndex = zeros(1,nq);
    for i=1:2*L
        for j=1:2*L
            if done_flag((i-1)*(2*L)+j) == 0
                % locate the conjugate of yi and yj in the previous W
                ic = cIndex(i);
                jc = cIndex(j);
                Q((i-1)*(2*L)+j) = 0.5*(W(i)*W(j)+W(ic)*W(jc));
                done_flag((i-1)*(2*L)+j) = 1;
                % calculate the conjugate of element (i-1)*(2*L)+j in the
                % current Q
                ij_c = (ic-1)*(2*L)+j;
                conjIndex((i-1)*(2*L)+j) = ij_c;
                % w(ic,jc) = w(i,j)
                Q((ic-1)*(2*L)+jc) = Q((i-1)*(2*L)+j);
                icjc_c = (i-1)*(2*L)+jc;
                conjIndex((ic-1)*(2*L)+jc) = icjc_c;
                done_flag((ic-1)*(2*L)+jc) = 1;
                % calculate w(i,jc)
                Q((i-1)*(2*L)+jc) = 0.5*(W(i)*W(jc)+W(ic)*W(j));
                done_flag((i-1)*(2*L)+jc) = 1;
                ijc_c = (ic-1)*(2*L)+jc;
                conjIndex((i-1)*(2*L)+jc) = ijc_c;
                % w(ic,j) = w(i,jc)
                Q((ic-1)*(2*L)+j) = Q((i-1)*(2*L)+jc);
                done_flag((ic-1)*(2*L)+j) = 1;
                icj_c = (i-1)*(2*L)+j;
                conjIndex((ic-1)*(2*L)+j) = icj_c;
            end
        end
    end
else % type 1 transformation, Q's size is 2*|W|^2
    nq = 2*(2*L)^2;
    done_flag = zeros(1,nq);
    Q = zeros(1,nq);
    conjIndex = zeros(1,nq);
    for i=1:2*L
        for j=1:2*L
            for k=1:2
                if(done_flag((i-1)*2*L+(j-1)*2*L+k) == 0)
                    u = k-1;
                    % calculate W(yi,yj,u|0) = 1/2*W(yi|u)*W(yj|0)
                    if u==0
                        idx_iu = i;
                    else
                        idx_iu = cIndex(i);
                    end
                    Q((i-1)*2*L+(j-1)*2*L+k) = 0.5*W(idx_iu)*W(j);
                    done_flag((i-1)*(2*L)+(j-1)*2*L+k) = 1;
                     % calculate the conjugate of element (ic-1)*(2*L)+(jc-1)*2*L+k in the current Q
                    ic = cIndex(i);
                    jc = cIndex(j);
                    iju_c = (ic-1)*(2*L)+(jc-1)*2*L+k;
                    conjIndex((i-1)*(2*L)+(j-1)*2*L+k) = iju_c;
                    % W(ic,j,uc) = W(i,j,u)
                    uc = mod(u+1,2); % in 0/1
                    kc = uc + 1; % convert to Matlab index
                    Q((ic-1)*2*L+(j-1)*2*L+kc) = Q((i-1)*2*L+(j-1)*2*L+k);
                    done_flag((ic-1)*(2*L)+(j-1)*2*L+kc) = 1;
                    icjuc_c = (i-1)*2*L+(jc-1)*2*L+kc;
                    conjIndex((ic-1)*(2*L)+(j-1)*2*L+kc) = icjuc_c;
                    % calculate W(ic,jc,u|0)
                    if u==0
                        idx_icu = ic;
                    else
                        idx_icu = cIndex(ic);
                    end
                    Q((ic-1)*2*L+(jc-1)*2*L+k) = 0.5*W(idx_icu)*W(jc);
                    done_flag((ic-1)*(2*L)+(jc-1)*2*L+k) = 1;
                    icjcu_c = (i-1)*(2*L)+(j-1)*2*L+k;
                    conjIndex((ic-1)*(2*L)+(jc-1)*2*L+k) = icjcu_c;
                    % W(i,jc,uc|0) = W(ic,jc,u|0)
                    Q((i-1)*2*L+(jc-1)*2*L+kc) =  Q((ic-1)*2*L+(jc-1)*2*L+k);
                    done_flag((i-1)*(2*L)+(jc-1)*2*L+kc) = 1;
                    ijcuc_c = (ic-1)*2*L+(j-1)*2*L+kc;
                    conjIndex((i-1)*(2*L)+(jc-1)*2*L+kc) = ijcuc_c;
                    % calculate W(yi,yj,kc|0)
                    if uc==0
                        idx_iuc = i;
                    else
                        idx_iuc = cIndex(i);
                    end
                    Q((i-1)*2*L+(j-1)*2*L+kc) = 0.5*W(idx_iuc)*W(j);
                    done_flag((i-1)*(2*L)+(j-1)*2*L+kc) = 1;
                    ijuc_c = (ic-1)*2*L+(jc-1)*2*L+kc;
                    conjIndex((i-1)*2*L+(j-1)*2*L+kc) = ijuc_c;
                    % W(ic,j,u | 0) = W(i,j,uc|0);
                    Q((ic-1)*2*L+(j-1)*2*L+k) = Q((i-1)*2*L+(j-1)*2*L+kc)
                    done_flag((ic-1)*2*L+(j-1)*2*L+k) = 1;
                    icju_c = (i-1)*2*L+(jc-1)*2*L+k;
                    conjIndex((ic-1)*2*L+(j-1)*2*L+k) = icju_c;
                    % calculate W(ic,jc,uc | 0)
                    if uc==0
                        idx_icuc = ic;
                    else
                        idx_icuc = cIndex(ic);
                    end
                    Q((ic-1)*2*L+(jc-1)*2*L+kc) = 0.5*W(idx_icuc)*W(jc);
                    done_flag((ic-1)*2*L+(jc-1)*2*L+kc) = 1;
                    icjcuc_c = (i-1)*2*L+(j-1)*2*L+kc;
                    conjIndex((ic-1)*2*L+(jc-1)*2*L+kc) = icjcuc_c;
                    % W(i,jc,u | 0) = W(ic,jc,uc | 0)
                    Q((i-1)*2*L+(jc-1)*2*L+k) = Q((ic-1)*2*L+(jc-1)*2*L+kc);
                    done_flag((i-1)*2*L+(jc-1)*2*L+k) = 1;
                    ijcu_c = (ic-1)*2*L+(j-1)*2*L+k;
                    conjIndex((i-1)*2*L+(jc-1)*2*L+k) = ijcu_c;
                end
            end
        end
    end
end


end

