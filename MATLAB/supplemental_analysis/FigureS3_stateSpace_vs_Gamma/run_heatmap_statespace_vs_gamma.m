%% =================================== Figure S3 ========================================================
% This file generates a state-space starting position vs gamma (CM activation effect) heatmap (Figure S3)
% =======================================================================================================

%% Critical points from data (see mouse_load_data from main code)
c1 = 1; c3 = 0; c2_mRNA = 0.6526; c2_miRNA = 0.7144;

%% Range of gamma and starting position values:
num_sim_runs = 100;          % Since this is a stochastic simulation requires multiple runs
gamma_values = 0.24:-0.01:0; % CM activation parameter range [0,0.24)

% Initial condition-starting point of trajectories to test
x0_mRNA = linspace(c3, c1, 9);
x0_miRNA = linspace(c3, c1, 9);
x0_mRNA = sort([x0_mRNA, c2_mRNA,c2_miRNA]);
x0_miRNA = sort([x0_miRNA, c2_mRNA,c2_miRNA]);

%% Store info
num_starting_pos = length(x0_mRNA);
num_gamma_values = length(gamma_values);

tot_num_cured_mRNA = zeros(num_starting_pos,num_gamma_values);
tot_num_cured_miRNA = zeros(num_starting_pos,num_gamma_values);
tot_num_cured_both = zeros(num_starting_pos,num_gamma_values);

%% Generate the heatmap
for i=1:num_sim_runs

    stateSpace_vs_diseaseGamma;

    %update counter:
    tot_num_cured_mRNA = tot_num_cured_mRNA + counter_cured_mRNA;
    tot_num_cured_miRNA = tot_num_cured_miRNA + counter_cured_miRNA;
    tot_num_cured_both = tot_num_cured_both + counter_cured_both;

    % if 90% or more are cured, then store a 1 to represent a health state
    cured_mRNA_final = double((tot_num_cured_mRNA./num_sim_runs) >= 0.90);
    cured_miRNA_final = double((tot_num_cured_miRNA./num_sim_runs) >= 0.90);    
    cured_both_final = double((tot_num_cured_both./num_sim_runs) >= 0.90);
end

%% Heatmaps for x0 vs gamma: Heatmap colors are defined as follows 
%    * Black for disease state = 0 
%    * Green for health state = 1

%===============================================
% mRNA (Figure S3 top panel)
%===============================================
figure;
h = heatmap(x0_mRNA, gamma_values, cured_mRNA_final.');
h.Colormap = [0 0 0; 
        0 0.7 0]; 
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';

%find the indices that correspond to c1, c2,c3:
[d,c2_mRNA_idx] = min(abs(x0_mRNA-c2_mRNA));
[d2,c2_miRNA_idx] = min(abs(x0_miRNA-c2_miRNA));
critical_pts_indx = [1,c2_mRNA_idx,c2_miRNA_idx,length(x0_mRNA)];

%Custom label the critical points
custom_labels = repmat({''}, length(x0_mRNA), 1);  
custom_labels{critical_pts_indx(1)} = 'c_3';
custom_labels{critical_pts_indx(2)} = 'c_{2,mRNA}';  
custom_labels{critical_pts_indx(3)} = 'c_{2,miR}';
custom_labels{critical_pts_indx(4)} = 'c_1'; 
h.XDisplayLabels = custom_labels;

xlabel('State space starting position')
ylabel('\gamma')
title('mRNA: x_0 vs \gamma')

%===============================================
% miRNA (Figure S3 bottom panel)
%===============================================
figure;
h = heatmap(x0_miRNA, gamma_values, cured_miRNA_final.');
h.Colormap = [0 0 0; 
        0 0.7 0]; 
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';

%find the indices that correspond to c1, c2,c3:
[d,c2_mRNA_idx] = min(abs(x0_mRNA-c2_mRNA));
[d2,c2_miRNA_idx] = min(abs(x0_miRNA-c2_miRNA));
critical_pts_indx = [1,c2_mRNA_idx,c2_miRNA_idx,length(x0_mRNA)];

%Custom label the critical points
custom_labels = repmat({''}, length(x0_mRNA), 1);  % Start with all labels blank
custom_labels{critical_pts_indx(1)} = 'c_3';
custom_labels{critical_pts_indx(2)} = 'c_{2,mRNA}';  
custom_labels{critical_pts_indx(3)} = 'c_{2,miR}';
custom_labels{critical_pts_indx(4)} = 'c_1'; 
h.XDisplayLabels = custom_labels;

xlabel('State space starting position')
ylabel('\gamma')
title('miRNA: x_0 vs \gamma')