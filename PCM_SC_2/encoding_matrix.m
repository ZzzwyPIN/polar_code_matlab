%%%%% LL Feb. 19, 2013
%%%% Encoding generating matrix GN
%%%% GN = BN * F@n
% n is the level of the encoding ==> block length N=2^n
function GN = encoding_matrix(n)

    % block length of the polar code
    N = 2^n;    
    F = [1 0; 1 1];

    % obtain BN through recursive relation BN = RN*(I2 @ BN/2)
    Bpre = 1;
    Fpre = 1;
    for i=1:n
        Ni = 2^i;
        B = zeros(Ni, Ni);
        Fn = zeros(Ni, Ni);
        % reverse shuffle matrix R
        R = zeros(Ni, Ni);
        for j=1:Ni/2
            R(2*j-1,j) = 1;
            R(2*j, Ni/2+j) = 1;
        end
        % recursive B
        B = R*kron(eye(2),Bpre);
        Bpre = B;
        % recursive FN
        Fn = kron(F, Fpre);
        Fpre = Fn;
    end
    
    GN = B*Fn;
end