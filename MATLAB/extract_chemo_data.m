function out = extract_chemo_data(data, ids, id_col, time_col, value_col, Yn_scaling, time_window, max_pts)

    out = extract_data(data, ids, id_col, time_col, value_col, Yn_scaling, time_window, max_pts);

    num_mice = length(ids);

    out.W0 = NaN(1,num_mice);               % start of treatment(W0)
    out.W0_idx = NaN(1,num_mice);
    out.max_pc_val = NaN(1,num_mice);       % max/peak response
    out.max_idx = NaN(1,num_mice);
    out.postchemo = NaN(max_pts,num_mice);  % post-chemo data points
    out.postchemo_pc_val = [];

    for i = 1:num_mice
        t = out.time(:,i); t = t(~isnan(t));
        v = out.value(:,i); v = v(~isnan(v));

        % --- Start of treatment (W0) ---
        idx0 = find(t==0,1);
        if ~isempty(idx0)
            out.W0(i) = v(idx0);
            out.W0_idx(i) = idx0;
        end

        % --- Post-chemo ---
        mask = t >= 0;        %post-chemo index
        pc_vals = v(mask);    %post-chemo pc-values

        if ~isempty(pc_vals)

            % find the max/peak value for each mouse trajectory
            [mx, idx] = max(pc_vals);
            orig_idx = find(mask);

            out.max_pc_val(i) = mx;
            out.max_idx(i) = orig_idx(idx);
            out.postchemo(1:length(pc_vals),i) = pc_vals;

            out.postchemo_pc_val = [out.postchemo_pc_val; pc_vals];
        end
    end
end

