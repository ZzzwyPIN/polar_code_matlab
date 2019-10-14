%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%This function produce the pure generation matrix, which not contain the
%%%bit-reversal matrix Bn
function G = spc_encoding(n)
    % block length of the polar code  
    F = [1 0; 1 1];
    Fpre = 1;
    for i=1:n
        % recursive FN
        Fn = kron(F, Fpre);
        Fpre = Fn;
    end
    G = Fn;
end



