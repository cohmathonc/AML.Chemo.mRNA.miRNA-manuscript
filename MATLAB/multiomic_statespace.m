% This function generates the 2D multiomic state-space (Figure 2A).
% The multiomic state-space is constructed by pairing the mRNA and miRNA state-space value based on time point 
% for each mouse in the 2018 untreated CM mice cohort (black) and chemotherapy treated mice (magenta).

function multiomic_statespace(data,params)

    % ---- Extract post-chemotherapy data ----
    mRNA_postchemo = data.mRNA.postchemo(:);
    miRNA_postchemo = data.miRNA.postchemo(:);

    mask_post = ~isnan(mRNA_postchemo) & ~isnan(miRNA_postchemo);
    mRNA_postchemo = mRNA_postchemo(mask_post);
    miRNA_postchemo = miRNA_postchemo(mask_post);

    % ---- Extract 2018 CM data ----
    mRNA_2018 = data.AML2018.mRNA.value(:);
    miRNA_2018 = data.AML2018.miRNA.value(:);

    mask_aml = ~isnan(mRNA_2018) & ~isnan(miRNA_2018);
    mRNA_2018 = mRNA_2018(mask_aml);
    miRNA_2018 = miRNA_2018(mask_aml);

    % ---- Regression (post-chemo) ----
    mdl = fitlm(mRNA_postchemo, miRNA_postchemo);
    [r, p] = corr(mRNA_postchemo, miRNA_postchemo);
    r2 = mdl.Rsquared.Ordinary;

    % Extract the slope (coefficient of x)
    slope_fit_chemo = mdl.Coefficients.Estimate(2);
    
    % Get prediction and confidence intervals
    %x_fit = linspace(min(mRNA_postchemo), max(mRNA_postchemo), 100)';
    x_fit = linspace(min(mRNA_postchemo), 1.1, 100)';
    [y_fit, y_ci] = predict(mdl, x_fit);
       
    % ---- Regression (AML 2018) ----
    mdl_AML = fitlm(mRNA_2018, miRNA_2018);
    [r_AML, p_AML] = corr(mRNA_2018, miRNA_2018);
    r2_AML = mdl_AML.Rsquared.Ordinary;

    % Extract the slope (coefficient of x)
    slope_fit_2018 = mdl_AML.Coefficients.Estimate(2);
    
    % Get prediction and confidence intervals
    x_fit_AML = linspace(min(mRNA_2018), max(mRNA_2018), 100)';
    [y_fit_AML, y_ci_AML] = predict(mdl_AML, x_fit_AML);

    % ---- Print stats ----
    fprintf('Post-chemo r=%.4f, R^2=%.4f, p=%.4f\n',  r, r2, p);
    fprintf('Slope of post-chemo fit line: %.4f\n', slope_fit_chemo);
    fprintf('AML r=%.4f, R^2=%.4f, p=%.4f\n', r_AML, r2_AML, p_AML);
    fprintf('Slope of AML fit line: %.4f\n', slope_fit_2018);

    % ---- Plot ----
    fig = figure; 
    hold on; box on;
    scatter(mRNA_2018, miRNA_2018, 30, 'k', 'filled', 'DisplayName','Untreated');
    scatter(mRNA_postchemo, miRNA_postchemo, 30, 'm', 'filled', 'DisplayName','Post-Chemo');
    % fit lines
    plot(x_fit, y_fit, 'm-', 'LineWidth',2, 'HandleVisibility','off');
    plot(x_fit_AML, y_fit_AML, 'k-', 'LineWidth',2, 'HandleVisibility','off');
    % 95% confidence intervals
    fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], [1 0 1], 'FaceAlpha',0.1, 'EdgeColor','none','HandleVisibility','off');
    fill([x_fit_AML; flipud(x_fit_AML)], [y_ci_AML(:,1); flipud(y_ci_AML(:,2))], [0 0 0], 'FaceAlpha',0.1, 'EdgeColor','none','HandleVisibility','off');
    % Axis formatting
    c_mRNA = data.mRNA.c;
    c_miRNA = data.miRNA.c;
    xlim([-0.2 1.2])
    ylim([-0.3 1.5])
    set(gca,'xtick',[c_mRNA(3) c_mRNA(2) c_mRNA(1)],'XTickLabel',{'$c_3$','$c_2$','$c_1$'},'TickLabelInterpreter','latex')
    set(gca,'ytick',[c_miRNA(3) c_miRNA(2) c_miRNA(1)], 'YTickLabel',{'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'}, 'TickLabelInterpreter','latex')
    xlabel('mRNA')
    ylabel('miRNA')
    title('Multiomic state-space')
    legend('Location','northwest');

    % ---- Save figure ----
    print(fig,'multiomic_statespace_Fig2A.svg','-dsvg')
    exportgraphics(fig,'multiomic_statespace_Fig2A.png','Resolution',600,'BackgroundColor','none')

end

