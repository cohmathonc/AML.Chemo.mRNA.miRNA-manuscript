% This function computes the mean-squared displacement (MSD) for each mRNA or miRNA mouse trajectory (Figure S1)
% Computation: MSD for each mouse transcriptome trajectory (mRNA or miRNA) prior to treatment (from T0-W0) was calculated. Since our potential changes after treatment, we cannot do MSD after W0.
% Then a linear fit is applied to the MSD and the slope of the linear regression is computed.
% The average slope (positive slopes only) of the linear regressions is used to estimate the diffusion coefficient for mRNA or miRNA trajectories. 

function result = compute_msd(data, type, save_name)

    % Get time and pc_value data
    time_matrix = data.time;
    pc_value_matrix = data.value;
    ids = data.ids;
    keepIdx = data.keepIdx; 

    num_mice = size(time_matrix,2);

    % Storage
    max_pts = size(time_matrix,1);

    MSD_matrix = NaN(max_pts-1, num_mice);
    time_lag_matrix = NaN(max_pts-1, num_mice);
    num_time_lags = zeros(1, num_mice);
    valid_ids = [];

    j = 1;

    % ===== Extract data points from T0 → W0 =====
    for i = 1:num_mice

        t = time_matrix(:,i);
        v = pc_value_matrix(:,i);

        t = t(~isnan(t));
        v = v(~isnan(v));

        if isempty(t)
            continue
        end

        % Only use trajectories that include pre-treatment (time points before W0 = 0)
        if any(t < 0)

            % Identify W0
            idx_W0 = find(t == 0, 1);

            if isempty(idx_W0)
                continue
            end

            % Truncate T0 → W0
            t_trunc = t(1:idx_W0);
            v_trunc = v(1:idx_W0);

            % Compute lag
            time_lag = abs(t_trunc(2:end) - t_trunc(1));

            % Compute MSD
            MSD = zeros(length(v_trunc)-1,1);

            for k = 1:length(v_trunc)-1
                displacement_sqrd = (v_trunc(2:k+1) - v_trunc(1)).^2;
                MSD(k) = mean(displacement_sqrd);
            end

            % Store
            time_lag_matrix(1:length(time_lag), j) = time_lag;
            MSD_matrix(1:length(MSD), j) = MSD;
            num_time_lags(j) = length(time_lag);

            valid_ids = [valid_ids, ids(i)];

            j = j + 1;
        end
    end

    % ===== Plot =====
    %fig = figure; 

    fig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    set(fig,'Color','w');

    clf;
    hold on;
   
    % Colormap
    cmap_full = colormap(hsv(12));
    cmap_full = cmap_full(keepIdx,:);

    subplot(1,2,1); hold on; box on; 
    title('Measured from data');
    subplot(1,2,2); hold on; box on;
    title('Linear fit');

    slopes = NaN(1,length(valid_ids));

    for k = 1:length(valid_ids)

        time_lag = time_lag_matrix(1:num_time_lags(k),k);
        msd = MSD_matrix(1:num_time_lags(k),k);

        % Fit
        p = polyfit(time_lag, msd, 1);
        slope = p(1);
        slopes(k) = slope;

        msd_fit = polyval(p, time_lag);

        % map ID → original index → color
        [~, mid_idx] = ismember(valid_ids(k), ids);
        color = cmap_full(mid_idx,:);

        % Plot measured
        subplot(1,2,1)
        plot(time_lag, msd, 'Color', color, 'LineWidth',2,'DisplayName', string(valid_ids(k)));

        % Plot fit
        subplot(1,2,2)
        plot(time_lag, msd_fit, 'Color', color,'LineWidth',2,'DisplayName', string(valid_ids(k)));

    end

    subplot(1,2,1)
    xlabel('Time (weeks)'); 
    ylabel('MSD');

    subplot(1,2,2)
    xlabel('Time (weeks)')

    sgtitle(type + " mean-squared displacement (MSD) analysis of mice trajectories")
    legend('Location','eastoutside')


    %% ===== Estimate diffusion coefficient - avg slope =====
    slopes_div2 = slopes / 2; % we use this slope based on Eistein's relation

    % Only positive slopes (important!)
    pos_slopes = slopes(slopes > 0);
    pos_slopes_div2 = slopes_div2(slopes_div2 > 0);

    % Store results
    result.slopes = slopes;
    result.slopes_div2 = slopes_div2;

    result.avg_slope = mean(slopes, 'omitnan');
    result.diffusion = mean(slopes_div2, 'omitnan');

    result.avg_positive_slope = mean(pos_slopes, 'omitnan');
    result.diffusion_positive = mean(pos_slopes_div2, 'omitnan'); % use only positive slopes to define the diffusion 
    result.valid_ids = valid_ids;

    % ===== Save figure =====
    save(save_name + "_results.mat", 'result')

    print(fig, save_name + "_FigS1" + ".svg", '-dsvg')
    exportgraphics(fig, save_name + "_FigS1"+ ".png",'Resolution',600,'BackgroundColor','none')

end

