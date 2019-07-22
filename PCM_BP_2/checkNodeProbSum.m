% this function calculates  \prod\alpha_i * phi( \sum phi(\beta_i) )
function sum1 = checkNodeProbSum(lr,log_flag)
    g_bp_max = 30;
    g_bp_min=phiBP(g_bp_max);
    n_lr = length(lr);
    sign_lr = 1;
    sum1 = 0;
%     lr_log=0;
    for i=1:n_lr
        if (log_flag==0)
            if (lr(i)==0)
                lr_log=-g_bp_max;
            else
                lr_log=log(lr(i));
            end
        else
            lr_log = lr(i);
        end
        if(lr_log < -30)
           lr_log=-g_bp_max;
        end
        if(lr_log < 0)
            sign_lr = sign_lr * (-1);
            lr_log = abs(lr_log);
        end
        % if lr_log is 0, assign the smallest lr to it
        if (lr_log == 0)
            % lr_log is 0, the result should be 0
%             disp('lr_log is 0')
%             lr_log = g_bp_min;
            sum1 = 0;
            return;
        end
%         tmp=0;
        tmp = phiBP(lr_log);
        % if tmp is infinite, assign the largest lr to it
        if isnan(tmp)
%             disp('inside summation is Inf');
            tmp = 0;
        end
        if isinf(tmp)
%             disp('inside summation is Inf');
            tmp = g_bp_max;
        end
        sum1 = sum1 + tmp;
    end
        if (sum1 == 0)
            sum1 = g_bp_min;
        end
        if (isinf(sum1)) 
            sum1 = g_bp_max;
        end
        sum1 = phiBP(sum1);
        if isnan(sum1)
    %         disp('NaN is encountered');
            sum1 = 0;
        end
        if isinf(sum1)
    %         disp('Inf is encountered');
            sum1 = g_bp_max;
        end
    sum1 = sign_lr * sum1;
    return
end