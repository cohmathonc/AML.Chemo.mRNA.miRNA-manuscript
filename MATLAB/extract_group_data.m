function out = extract_group_data(data, group_name, id_col, time_col, value_col, Yn_scaling)

    ids = unique(data.(id_col)(contains(data.treat_group_x, group_name)));

    num_mice = length(ids);
    out.time = NaN(11, num_mice);
    out.value = NaN(11, num_mice);
    out.group_pc_val = [];

    for i = 1:num_mice
        ind = data.(id_col) == ids(i);

        t = round(data.(time_col)(ind));
        v = Yn_scaling(data.(value_col)(ind));
        
        % Sort time in increasing order
        [t, s] = sort(t);
        v = v(s);

        % For 2018 CM mice with mouseIDs 3336 and 3357 we only want t<8 and t<5 respectively as in previous publications (Rockne et al.Cancer Res. 2020 & Frankhouser et al. Sci Adv. 2022)
        if ids(i) == 3336
            mask = t < 8;
            t = t(mask); v = v(mask);
        elseif ids(i) == 3357
            mask = t < 5;
            t = t(mask); v = v(mask);
        end

        % Store data for all mice:
        out.time(1:length(t), i) = t;
        out.value(1:length(v), i) = v;

        out.group_pc_val = [out.group_pc_val; v];
    end
end

