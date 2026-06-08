function c_scaled = apply_scaling(raw_c, xmin, xmax)
    
    %rescale state-space to 0-1
    Yn = @(x)(x - xmin) ./ (xmax - xmin);
    c = Yn(raw_c);
    
    %rescale critial points to 0-1
    xmin_new = c(3);
    xmax_new = c(1);

    Yn_new = @(x)(x - xmin_new) ./ (xmax_new - xmin_new);

    c_scaled = Yn_new(c); %c3 = 0 and c=1
end

