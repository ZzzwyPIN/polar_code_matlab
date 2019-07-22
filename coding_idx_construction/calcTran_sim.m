function [Q] = calcTran_sim( W, type_flag )
% calcTran_sim calculate the transition probability of the one-step channel transformation
%   W is the original channel transition probabilities
% type_flag: 0, type 0 transform; 1: type 1 transform
    L = length(W)/2;
    if type_flag == 0 % type 0 transformation, Q's size squares that of W
            nq = (2*L)^2;
            Q = zeros(1,nq);
            for i=1:L
                for j=1:L
                    %  the conjugate of yi and yj in the previous W is i+L and j+L
                    Q((i-1)*(2*L)+j) = 0.5*(W(i)*W(j)+W(i+L)*W(j+L));
                    % calculate w(ic,j)
                    Q((i-1)*(2*L)+j+2*L*L) = 0.5*(W(i+L)*W(j)+W(i)*W(j+L));
                    % w(i,jc)=w(ic,j)
                    Q((i-1)*(2*L)+j+L) = Q((i-1)*(2*L)+j+2*L*L);
                    % w(ic,jc) = w(i,j)
                    Q((i-1)*(2*L)+j+L+2*L*L) = Q((i-1)*(2*L)+j);
                end
            end
    else % type 1 transformation, Q's size is 2*|W|^2
        nq = 2*(2*L)^2;
        Q = zeros(1,nq);
        % W(yi,yj,u1 | 0 ) = 1/2 * W(yi|u1) * W(yj|0)
        for i=1:L
            for j=1:L
                    % calculate W(yi,yj,0 | 0)
                    Q((i-1)*4*L+(j-1)*2+1) = 0.5*W(i)*W(j);
                    % calculate W(yi,yj,1 | 0)
                    Q((i-1)*4*L+(j-1)*2+2) = 0.5*W(i+L)*W(j);
                    % calculate W(yi,yjc,0 | 0)
                    Q((i-1)*4*L+(j-1)*2+2*L+1) = 0.5*W(i)*W(j+L);
                    % calculate W(yi,yjc,1 | 0)
                    Q((i-1)*4*L+(j-1)*2+2*L+2) = 0.5*W(i+L)*W(j+L);
                    
                    % W(yic,yjc,0) = W(yi,yjc,1|0)
                    Q((i-1)*4*L+(j-1)*2+1+4*L*L) = Q((i-1)*4*L+(j-1)*2+2*L+2);
                    % W(yic,yjc,1 | 0) = W(yi,yjc,0 | 0)
                    Q((i-1)*4*L+(j-1)*2+2+4*L*L) = Q((i-1)*4*L+(j-1)*2+2*L+1);
                    % W(yic,yj,0 | 0) = W(yi,yj,1 | 0)
                    Q((i-1)*4*L+(j-1)*2+2*L+1+4*L*L) = Q((i-1)*4*L+(j-1)*2+2);
                    % W(yic,yj,1 | 0) = W(yi,yj,0 | 0)
                    Q((i-1)*4*L+(j-1)*2+2*L+2+4*L*L) = Q((i-1)*4*L+(j-1)*2+1);

                    % the corresponding conjugate symbol should be directly
                    % 4*L*L away
%                     % W(yic,yjc,0) = W(yi,yjc,1|0)
%                     Q((i-1)*4*L+(j-1)*2+2*L+1+4*L*L) = Q((i-1)*4*L+(j-1)*2+2*L+2);
%                     % W(yic,yjc,1 | 0) = W(yi,yjc,0 | 0)
%                     Q((i-1)*4*L+(j-1)*2+2*L+2+4*L*L) = Q((i-1)*4*L+(j-1)*2+2*L+1);
%                     % W(yic,yj,0 | 0) = W(yi,yj,1 | 0)
%                     Q((i-1)*4*L+(j-1)*2+1+4*L*L) = Q((i-1)*4*L+(j-1)*2+2);
%                     % W(yic,yj,1 | 0) = W(yi,yj,0 | 0)
%                     Q((i-1)*4*L+(j-1)*2+2+4*L*L) = Q((i-1)*4*L+(j-1)*2+1);
             end
        end
    end
end


