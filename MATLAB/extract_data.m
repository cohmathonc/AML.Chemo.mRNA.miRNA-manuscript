function out = extract_data(data, ids, id_col, time_col, value_col, transform, time_window, max_pts)

    num_mice = length(ids);
    out.time = NaN(max_pts, num_mice);
    out.value = NaN(max_pts, num_mice);
    out.time_array = [];
    out.value_array = [];

    for i = 1:num_mice
        ind = data.(id_col) == ids(i);

        t = round(data.(time_col)(ind));         %time
        v = transform(data.(value_col)(ind));    %value (ex: ckit, pc value)
        
        % Sort time in increasing order
        [t, s] = sort(t);
        v = v(s);
        
        % only want time between -3 weeks and 10 weeks
        mask = (t >= time_window(1)) & (t <= time_window(2));
        t = t(mask);
        v = v(mask);
        
        % Store data for all mice:
        out.time(1:length(t), i) = t;
        out.value(1:length(v), i) = v;

        out.time_array = [out.time_array; t];
        out.value_array = [out.value_array; v];
    end
end