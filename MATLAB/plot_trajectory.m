function [t_ave, y_ave] = plot_trajectory(data, color_avg, ylab, figure_title, figure_name, params, type, varargin)

    time_matrix = data.time;
    pc_value_matrix = data.value;
    t_array = data.time_array;
    pc_val_array = data.value_array;

    num_mice = size(time_matrix,2);
    num_pts = sum(~isnan(time_matrix),1);

    fig = figure; hold on;
    box on;

    % Individual trajectories
    for i = 1:num_mice
        t = time_matrix(1:num_pts(i), i);
        v = pc_value_matrix(1:num_pts(i), i);

        plot(t, v, 'Color',[.75 .75 .75], 'LineWidth',1.2);
        scatter(t, v, 'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','none');
    end

    % ===== Average trajectory =====
    t_ave = -3:params.tf;
    y_ave = NaN(size(t_ave));

    for k = 1:length(t_ave)
        mask = (t_array == t_ave(k));
        if any(mask)
            y_ave(k) = mean(pc_val_array(mask));
        end
    end

    valid = ~isnan(y_ave);
    t_ave = t_ave(valid);
    y_ave = y_ave(valid);

    plot(t_ave, y_ave, 'Color', color_avg, 'LineWidth',2);
    scatter(t_ave, y_ave, 'MarkerFaceColor', color_avg, 'MarkerEdgeColor', color_avg);

    % Save averages
    save(figure_name + "_avg.mat", 't_ave', 'y_ave');

    % ===== Critical points (optional) =====
    if ~isempty(varargin)
        c = varargin{1};

        yline(c(1),':k','LineWidth',2);
        yline(c(2),':','Color',[0.3 0.3 0.3],'LineWidth',2);
        yline(c(3),':r','LineWidth',2);

        ylim([c(3)-0.2 c(1)+0.1]);

        if strcmp(type, 'mRNA')
            labels = {'$c_3$','$c_2$','$c_1$'};
        elseif strcmp(type, 'miRNA')
            labels = {'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'};
        else
            labels = {'c_3','c_2','c_1'};
        end

        set(gca,'ytick',[c(3) c(2) c(1)],'YTickLabel',labels,'TickLabelInterpreter','latex')
    end

    % ===== Treatment shading =====
    yl = ylim;   % get [ymin ymax]
    fill([0 1 1 0], [yl(1) yl(1) yl(2) yl(2)],'m', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    % ===== Labels =====
    xlabel('Time (weeks)')
    ylabel(ylab)
    title(figure_title)

    xlim([-3.5 10.5])
    xticks(-3:1:10)

    % Special formatting for cKit (0–100 axis)
    if strcmp(ylab,'cKit%')
        yticks(0:10:100)
        labels = string(0:10:100);
        labels(mod(0:10:100,20)~=0) = '';
        yticklabels(labels)
    end

    % ===== Save figure =====
    print(fig, figure_name + ".svg", '-dsvg')
    exportgraphics(fig, figure_name + ".png",'Resolution',600,'BackgroundColor','none')

end

