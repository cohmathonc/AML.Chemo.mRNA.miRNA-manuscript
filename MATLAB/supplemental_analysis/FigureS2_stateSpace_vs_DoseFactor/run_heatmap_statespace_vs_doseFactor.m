%% =================================== Figure S2 =====================================================
% This file generates the state-space starting position vs dose factor (alpha) heatmap (Figure S2)
% Note: This file requires the user input
% *T (final time): for 10 weeks post-Rx use T = 13 and for 20 weeks post-Rx use T = 23
% *gamma (CM activation parameter): parameter range is [0,0.25). Pick a value.
% ===================================================================================================

%% Critical points come from data (see mouse_load_data.m from main code)
c1 = 1; c3 = 0; c2_mRNA = 0.6526; c2_miRNA = 0.7144;

%% User input:

num_sim_runs = 100;    %since this is a stochastic simulation requires multiple runs
T = 13;
gamma = 0.07;
dose_factors = linspace(3,1,11); % this will test 11 dose factors from 1 to 3 

% INPUT-Initial condition-starting point of trajectories to test
x0_mRNA = linspace(c3, c1, 9);
x0_miRNA = linspace(c3, c1, 9);
x0_mRNA = sort([x0_mRNA, c2_mRNA,c2_miRNA]);
x0_miRNA = sort([x0_miRNA, c2_mRNA,c2_miRNA]);

%% Store info
num_starting_pos = length(x0_mRNA);
num_dose_factors = length(dose_factors);

cured_both_counter = zeros(num_starting_pos,num_dose_factors);
cured_both_percent = zeros(num_starting_pos,num_dose_factors);

%% Generate the heatmap
for i=1:num_sim_runs
    stateSpace_vs_dose; % preRx data beta from T0-W0 and postRx data beta from W0-10
  
    %update counter:
    cured_both_counter = cured_both_counter + cured_both;
    cured_both_percent = cured_both_counter./num_sim_runs;

    %if 90% or more, then store a 1
    cured_both_final = double(cured_both_percent >= 0.90);
end

%% Heatmap for x0 vs alpha (dose factors):

figure;
h = heatmap(x0_mRNA, dose_factors, cured_both_final.');
h.Colormap = [0 0 0; %Black for disease = 0
        0 0.7 0]; %Green for health = 1
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
ylabel('\alpha')
title(sprintf('Both Trajectories (\\gamma = %.2f)', gamma))

%% Plot E for different dose factors
figure;
plot(stime,E_matrix, 'LineWidth',1); 
grid on;
xlabel('Time');
ylabel('E(t)')
title('Effect per dose factor \gamma_{d}');
