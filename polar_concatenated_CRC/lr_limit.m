function lr_r = lr_limit( lr, lr_1,lr_2, z_type)
% This function limit the lr value from possible numerical overflow
% lr : the result lr value
% lr_1, lr_2: the lr value used to calculate lr
% z_type: lr is the upper or lower at the Z-connection.0: upper, 1:lower

r_max = realmax; % find the largest possible value of Matlab
r_min = realmin; % find the smallest possible value of Matlab

lr_r = lr;

if lr == 0
    lr_r = r_min;
end

if isinf(lr)
    lr_r = r_max;
end

if isnan(lr)
    % upper
    if z_type == 0
        if(lr_1 < lr_2)
            lr_r = lr_1;
        else
            lr_r = lr_2;
        end
    end
    % if lower, no need to process lr.
    % the way lr is processed makes sure that there is no NaN for lower LR
end

end

