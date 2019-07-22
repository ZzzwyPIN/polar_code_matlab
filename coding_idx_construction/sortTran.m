function [Q,conjIndex] = sortTran( W,cIndex )
% sortTran sort the input W according to the LR values
%  The output Q has entries of y1,y2,..., yL,y1c,y2c,...,yLc with LRs
%  in the accending order

    nw = length(W);
    L = nw/2;
    lr = zeros(1,nw);
    done_flag = zeros(1,nw);
    for i=1:nw
        if(done_flag(i) == 0)
            ic = cIndex(i);
            % possible hardware improvement
            % to make the final yi have LR(yi) > 1
            % otherwise exchange yi and yic
            if(W(i) > W(ic))
                lr(i) = W(i) / W(ic);
            else
                lr(i) = W(ic) / W(i);
                tmp = W(i);
                W(i) = W(ic);
                W(ic) = tmp;
            end
            lr(ic) = 1/lr(i);
            done_flag(i) = 1;
            done_flag(ic) = 1;
        end
    end
    for i=1:nw
    end
    % sort W according to LR's 
end

