function [Q] = sortTran_sim( W)
% sortTran sort the input W according to the LR values
%  The output Q has entries of y1,y2,..., yL,y1c,y2c,...,yLc with LRs
%  in the accending order
    nw = length(W);
    L = nw/2;
    lr = zeros(1,L);
    lr = W(1:L) ./ W(L+1:nw);
    % select from each pair (y,yc) so that LR(y) > 1
    Is = find(lr < 1);
    wtmp = W(Is);
    W(Is) = W(Is+L);
    W(Is+L) = wtmp;
    lr(Is) = 1./lr(Is);
    [lrs,I] = sort(lr);
    Q = [W(I) W(I+L)];
end

