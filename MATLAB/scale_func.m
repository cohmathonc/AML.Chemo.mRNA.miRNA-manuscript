function scale_fn = scale_func(xmin, xmax, raw_c)

    scale_fn = @(x)(x - xmin) ./ (xmax - xmin);

end

