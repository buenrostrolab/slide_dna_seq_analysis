function[local_mins] = local_min_2d(x)

    local_mins = zeros(size(x));
    
    for i=2:size(x,1)-1;
        local_mins(i) = x(i) < x(i-1) & x(i) < x(i+1);
    end
end