function [Ixy, P] = calcMutual( Wyx )
% calculate the mutual information of Y and X and returns the probability
% of error
%   Wyx is the given transition probability of y given x
%   This function assumes a fixed pattern of Wyx: y1,y2,...,yL,
%   y1c,y2c,...,yLc where y1c is the conjugate symbol of y1

    nq = length(Wyx);
    P = sum(Wyx(nq/2+1:nq))/2;
    % calculate the mutual infomation
    % w(y)
    Wy = zeros(1,nq/2);
    Wy = (Wyx(1:nq/2) + Wyx(nq/2+1:nq))/2;
    % Hy
    Hy = -2*sum(Wy .* log2(Wy));
    % H(Y|X)
    Hyx = -sum(Wyx.*log2(Wyx));
    Ixy = Hy-Hyx;

end

