% This code runs the two-dimensional (2D) state-transition model 
% to study the temporal dynamics of peripheral blood mRNA and microRNA miRNA transcriptomes in a mouse model of AML
% Executing this code reproduces the analysis and figures generated for the manuscript

addpath('supplemental_analysis')

%% ====== Parameters ======

% Time:
params.t0 = 0; %initial time point (weeks)
params.tf = 20; %max final time point of all data (weeks)

% mRNA critical points and scaling based on state-space dimensions
params.mRNA.raw_cp = [0.5236, -84.3842,-243.8802]; %c1, c2, c3
params.mRNA.xmin = -205.9283; 
params.mRNA.xmax = 15.46172;

% miRNA critical points and scaling based on state-space dimensions
params.miRNA.raw_cp = [3.151,-13,-53.395]; %c1, c2, c3
params.miRNA.xmin = -29.15299; 
params.miRNA.xmax = 2.139824;

%% ====== Load the data ======
% Data is imported from the data file "m_miRNA_combined_all.csv" created from the analysis of the mRNA-seq and miRNA-seq data.
% See the "Chemo_AML_Analysis.R" file in the R folder.
data = load_mouse_data(params);

%% ====== Data trajectories (Figure 1) ======

% ------ Plot cKit% over time (Figure 1C) ------
plot_trajectory(data.ckit, [0.7 0 0],'cKit%','cKit%','cKit_Fig1C', params);

% ------ Transcriptome state-space trajectories (Figure 1D-E) ------
[t_ave_mRNA, y_ave_mRNA] = plot_trajectory(data.mRNA,[0.3059, 0.6549, 0.1804],'AML state-space','mRNA','mRNA_data_trajectories_Fig1D',params, 'mRNA', data.mRNA.c);
[t_ave_miRNA, y_ave_miRNA] = plot_trajectory(data.miRNA,[0.1294, 0.3725, 0.6039],'AML state-space','miRNA','miRNA_data_trajectories_Fig1E',params, 'miRNA', data.miRNA.c);

% ------ Average trajectories (Figure 1F) ------
fig = figure; 
ax = gca;
hold on; box on;
plot(t_ave_mRNA,y_ave_mRNA,'Color', [0.3059, 0.6549, 0.1804], 'LineWidth',2); hold on;
plot(t_ave_miRNA,y_ave_miRNA,'Color',[0.1294, 0.3725, 0.6039], 'LineWidth',2); hold on;
% Get mRNA critical points:
c = data.mRNA.c;    
% Define limits for treatment shaded region:
y_min = c(3)-0.2; 
y_max = c(1)+0.1;
ylim([y_min y_max])
xlim([-3.5 10.5])
fill([0 1 1 0], [y_min y_min y_max y_max],'m', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% Remove ticks:
xlabel('Time (weeks)')
xticks(-3:1:10)
ax.XAxisLocation = 'bottom';
set(ax, 'YTickLabel', []);
set(ax, 'YTick', []);
print(fig,'avg_data_trajectories_Fig1F.svg','-dsvg')
exportgraphics(fig,'avg_data_trajectories_Fig1F.png','Resolution',600,'BackgroundColor','none')

%% ====== Multiomic state-space and 2D state-transition model (Figure 2) ======

% Multiomic state-space with 2018 untreated CM mice (black) and chemotherapy treated mice (magenta) (Figure 2A):
multiomic_statespace(data, params);

% Perform an MSD analysis to identify diffusion coefficients for Langevin eqn in 2D model (Figure S1)
% MSD function is in "supplemental_analysis" folder
MSD_mRNA = compute_msd(data.mRNA, 'mRNA', 'MSD_mRNA');
MSD_miRNA = compute_msd(data.miRNA, 'miRNA', 'MSD_miRNA');

Diff.mRNA = MSD_mRNA.diffusion_positive;
Diff.miRNA = MSD_miRNA.diffusion_positive;

% 2D state-transition model (Figure 2B-E):
run_2D_stateTransition_model;

% Compute the treatment vector (Figure 2F):
Treatment_vector;