function GN = encoding_matrix2(n)

    % block length of the polar code
    N = 2^n;    
    %F = [1 0; 1 1];

    % obtain BN through recursive relation BN = RN*(I2 @ BN/2)
    Gpre = 1;
    %Fpre = 1;
    for i=1:n
        Ni = 2^i;
        G = zeros(Ni, Ni);
        %Fn = zeros(Ni, Ni);
        % reverse shuffle matrix R
        A = zeros(Ni, Ni);
        for j=1:Ni/2
            A(2*j-1,j) = 1;
            A(2*j,j) = 1;
            A(2*j, Ni/2+j) = 1;
        end
        % recursive B
        G = A*kron(eye(2),Gpre);
        Gpre = G;
        % recursive FN
        %Fn = kron(F, Fpre);
        %Fpre = Fn;
    end
    
    GN = G;
end