function [Q] = degrading_merge( W,u)
% merges the output of the channel W with a fixed output alphabet size u=2v
% the alphabet size of W is 2L

debugging = 0;

L = length(W)/2;
v = u/2;
if v >= L
    Q = W;
else
    deltaI = zeros(1,L-1);
    for i=1:L-1
        deltaI(i) = calcDeltaI(W(i),W(i+L),W(i+1),W(i+1+L));
    end
    while(L > v)
        % debugging
        if debugging
            if L==5
                L
            end
        end
        deltaI = deltaI(1:(L-1));
        min_idx_tmp = find(deltaI == min(deltaI));
        % if more than one min is found
        min_idx = min_idx_tmp(1);
        ap = W(min_idx) + W(min_idx+1);
        bp = W(min_idx+L) + W(min_idx+1+L);
        % symbols with indices smaller than min_idx stay at the same memory
        % elements
        Wtmp = zeros(1,2*(L-1));
        if(min_idx-1 > 0)
            deltaI(min_idx-1) = calcDeltaI(W(min_idx-1),W(min_idx-1+L),ap,bp);
            Wtmp(1:min_idx-1) = W(1:min_idx-1);
            Wtmp((1:min_idx-1)+L-1) = W((1:min_idx-1)+L);
        end
        % symbols with indices larger than min_idx move one position to the
        % left
        if(min_idx+2<=L)
            deltaI(min_idx) = calcDeltaI(ap,bp,W(min_idx+2),W(min_idx+2+L));
            Wtmp(min_idx+1:L-1) = W(min_idx+2:L);
            Wtmp((min_idx+1:L-1)+L-1) = W((min_idx+2:L)+L);
        end
        % merge y(min_idx) and y(min_idx+1) and replace these two symbols
        % by a new symbol z
        Wtmp(min_idx) = ap;
        Wtmp(min_idx+L-1) = bp;
        W(1:2*(L-1)) = Wtmp;
        % debugging
        if debugging
            checkSum = sum(W(1:2*(L-1))-Wtmp)
            if checkSum ~= 0
                W(1:2*(L-1))
                Wtmp
            end
        end
        deltaI(min_idx+1:L-2) = deltaI(min_idx+2:L-1);
        L = L-1;
    end
    Q=W(1:2*v);
end

end

