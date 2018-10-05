function [dec, dec_list] = polar_decoder_scl(llr, K, L)

N = length(llr);
n = log2(N);
[~, info_bit_idx, frozen_bit_flag] = polar_seq_gen(N, K);

% Initialization the decoding list
list_size = 1;
dec_list  = cell(1,L);
for listIdx = 1:list_size
    % Each list contains a LLR array, partial sum array and path metric(PM)
    dec_list{listIdx}.LLR = [NaN(N,log2(N)), llr(:)];
    dec_list{listIdx}.par_sum = NaN(N,log2(N)+1);
    dec_list{listIdx}.PM  = 0;
end

% pre-calculate code trellis
trellis = zeros(n, 2^(n-1));
for s = 1:n
    idx = repmat((0:2^(n-s)-1)*2^s,2^(s-1),1);
    trellis(s,:) = idx(:)'+ repmat(1:2^(s-1),1,2^(n-s));
end

% SCL decoding
for i = 1:N
    % update LLR before decoding bit i
    for listIdx = 1:list_size
%         fprintf('1,listIdx=%d\n', listIdx);
        LLR = dec_list{listIdx}.LLR;
        par_sum = dec_list{listIdx}.par_sum;
        % forward calculation of partial sum
        for s = 1:n
            idx = trellis(n-s+1,:);
            idx = idx(~(isnan(par_sum(idx,s)) | isnan(par_sum(idx+2^(n-s),s))) & isnan(par_sum(idx,s+1)));
            par_sum(idx,s+1) = mod(par_sum(idx,s)+par_sum(idx+2^(n-s),s), 2);
            par_sum(idx+2^(n-s),s+1) = par_sum(idx+2^(n-s),s);
        end
        % backward calculation of LLR
        for s = n:-1:1
            % perform f operation
            idx = trellis(n-s+1,:);
            idx = idx(~(isnan(LLR(idx,s+1)) | isnan(LLR(idx+2^(n-s),s+1))) & isnan(LLR(idx,s)));
            LLR(idx,s) = f(LLR(idx,s+1),LLR(idx+2^(n-s),s+1),0);
            % perform g operation
            idx = trellis(n-s+1,:);
            idx = idx(~isnan(par_sum(idx, s)) & isnan(LLR(idx+2^(n-s),s)));
            LLR(idx+2^(n-s),s) = g(par_sum(idx,s),LLR(idx,s+1),LLR(idx+2^(n-s),s+1));
            %fprintf(' LLR(idx+2^(n-s),s)=%d,idx+2^(n-s)=%d,par_sum(idx,s)=%d\n', LLR(idx+2^(n-s),s),idx+2^(n-s),par_sum(idx,s));
        end
        dec_list{listIdx}.LLR = LLR;
        dec_list{listIdx}.par_sum = par_sum;
%         fprintf('2,dec_list{listIdx}.LLR(2,1)=%d,listIdx=%d\n',dec_list{listIdx}.LLR(2,1), listIdx);
    end
    
    % update path metric and decoding list
   k =i-1;
	

k_rever = 0;
	

for m = 0:n-1
	

k_rever = bitor(k_rever,bitset(0,n-m,bitand(bitshift(1,m),k )));
	

end
	

i_0 = k_rever+1;

    if frozen_bit_flag(i)
        % frozen bit
        for listIdx = 1:list_size
            % ui 0
            if(0 ~= h(dec_list{listIdx}.LLR(i_0,1)))
                dec_list{listIdx}.PM = phi(dec_list{listIdx}.PM, dec_list{listIdx}.LLR(i_0,1), 0, 1);
            else % not update PM
            end
            dec_list{listIdx}.par_sum(i_0,1) = 0;
%             fprintf('dec_list{listIdx}.PM=%d\n', dec_list{listIdx}.PM);
        end
    else
        % info bit
        PM = zeros(1,2*list_size);
        for listIdx = 1:list_size
            PM(listIdx) = phi(dec_list{listIdx}.PM, dec_list{listIdx}.LLR(i_0,1), 0, 0);
            PM(listIdx+list_size) = phi(dec_list{listIdx}.PM, dec_list{listIdx}.LLR(i_0,1), 1, 0);
            dec_list{listIdx}.par_sum(i_0,1) = 0;
            dec_list{listIdx}.PM = PM(listIdx);
            dec_list{listIdx+list_size} = dec_list{listIdx};
            dec_list{listIdx+list_size}.par_sum(i_0,1) = 1;
            dec_list{listIdx+list_size}.PM = PM(listIdx+list_size);
        end
        if list_size >= L
%             fprintf('3,list_size=%d\n', list_size);
            [~,sort_idx] = sort(PM);
            dec_list = dec_list(sort_idx(1:L));
        else
            list_size = 2*list_size;
%             fprintf('3,list_size=%d\n', list_size);
        end
    end
end

dec = dec_list{1}.par_sum(info_bit_idx+1,1)';

end



function z = f(x,y,dec_type)
if dec_type == 0
    z = 2.*atanh(tanh(x/2).*tanh(y/2));
else
    z = sign(x).*sign(y).*min(abs(x),abs(y));
end
end

function z = g(u,x,y)
z = (1-2*u).*x + y;
end

function z = h(x)
z = 1/2*(1-sign(x));
end

function PM = phi(PM,LLR,u,dec_type)
if dec_type == 0
%     fprintf('PMini=%d,LLR=%d,u=%d\n',PM,LLR,u);
    if u ~= LLR
       PM = PM + log(1+exp(-(1-2*u)*LLR)); 
    else % PM = PM
    end
else
%     fprintf('PMini=%d,LLR=%d,u=%d\n',PM,LLR,u);
    PM = PM + ((1-2*u) ~= sign(LLR))*abs(LLR);
end

end