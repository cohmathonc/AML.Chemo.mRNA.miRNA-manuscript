%% ============================================================
%  2D STATE-TRANSITION MODEL (Figure 2)
%  Simulates stochastic dynamics in multiomic state-space
%  using a Langevin equation of motion
%% ============================================================

rng(1,'twister')
%rng('shuffle')

%% Simulation parameters

% ===== Time parameters =====
T = 13;
dt_sim = 0.1;
t = 0:dt_sim:T;
w0 = 3;      %start of treatment
wf = w0+1;   %end of treatment

% ===== Critical points =====
c1 = 1;
c3 = 0;
c2_mRNA  = data.mRNA.c(2);
c2_miRNA = data.miRNA.c(2);

% ===== Diffusion coefficient from MSD analysis =====
Dx_mRNA = MSD_mRNA.diffusion_positive;
Dx_miRNA = MSD_miRNA.diffusion_positive;
Diff = diag([sqrt(2*Dx_mRNA), sqrt(2*Dx_miRNA)]); % diffusion matrix

% ===== Defining the meshgrid ===== 
xmin = -1; 
xmax = 2; 
myx = linspace(xmin,xmax);
myy = myx;
[X,Y] = meshgrid(myx,myy);

% ===== Parameter that represents the effect of CM activation =====
gamma = 0.08;

%% Treatment parameters
alpha = 0.06;
alpha_treated = alpha;
alpha_untreated = 1; %alpha = 1 to get untreated conditions

% ===== Tensor operators psi and phi under treatment: cytarabine (50mg/kg/day; 5 days), daunorubicin (1.5mg/kg/day, 3 days) =====
psi = [1 1; 0 10];
phi = [5 0.1; 0 3];

% Under maximum treatment effect (E = Emax)
psi_maxEffect = psi;
phi_maxEffect = phi;

% WITHOUT treatment (identity matrix)
psi_untreated = [1 0; 0 1];
phi_untreated = [1 0; 0 1];

% ===== Chemo effect from PK-PD model =====
chemo_effect_PKPD; 
[effect_max, idx_max] = max(E); %maximum effect
rr = 1;                         %turn effect on/off

% ===== Initial condition (from data) =====
%Data starting points:
x0_mRNA = [0.6658, 0.4530, 0.6398, 0.6351, 0.5055, 0.5798];
x0_miRNA = [0.6740, 0.5148, 0.6541, 0.6442, 0.4393, 0.5853];

x0 = [x0_mRNA(2), x0_miRNA(2)]; % Initial condition-starting point of trajectories

%% Define untreated and treated potentials where SDE = [-dF/dx; -dF/dy]

% ===== Untreated potential & drift =====
F_untreated = (X.^2 + Y.^2).*(((X-c1).^2 + (Y-c1).^2) + gamma); %2D potential with wells at (0,0) and (1,1)
SDE = @(t, x) [-((2.*x(1).*(((x(1)-c1).^2 + (x(2)-c1).^2) + gamma)) + ((x(1).^2 + x(2).^2) .* (2.*(x(1) - c1)))); -((2.*x(2).*(((x(1)-c1).^2 + (x(2)-c1).^2)+ gamma)) + (2.*(x(2)-c1).*(x(1).^2 + x(2).^2)))];

% ===== Treated potential & drift =====
SDE_chemo = @(t,x) [-alpha.*(((2*psi(1,1).*x(1) + (psi(1,2) + psi(2,1)).*x(2)).*((phi(1,1).*(x(1)-c1).^2 + (phi(1,2) + phi(2,1)).*(x(1)-c1).*(x(2)-c1) + phi(2,2).*(x(2)-c1).^2) + gamma)) + (psi(1,1).*x(1).^2 + (psi(1,2) + psi(2,1)).*x(1).*x(2) + psi(2,2).*x(2).^2 + rr.*E(find((stime<=t),1,'last'))).*(2*phi(1,1).*(x(1)-c1) + (phi(1,2) + phi(2,1)).*(x(2)-c1))); -alpha.*((((psi(1,2) + psi(2,1)).*x(1) + 2*psi(2,2).*x(2)).*((phi(1,1).*(x(1)-c1).^2 + (phi(1,2) + phi(2,1)).*(x(1)-c1).*(x(2)-c1) + phi(2,2).*(x(2)-c1).^2) + gamma)) + ((psi(1,1).*x(1).^2 + (psi(1,2) + psi(2,1)).*x(1).*x(2) + psi(2,2).*x(2).^2 + rr.*E(find((stime<=t),1,'last'))).*((phi(1,2) + phi(2,1)).*(x(1)-c1) + 2*phi(2,2).*(x(2)-c1))))];

%% Simulate SDE

%  ===== Storage =====
%Store tensor elements
A11 = []; A12 = []; A21 = [];A22 = [];
B11 = []; B12 = []; B21 = [];B22 = [];
Alpha_array = [];
treatment_time = [];
effect_matrix = [];

x = zeros(2, length(t));  %store particle position
x(:, 1) = x0;             %initial position

% ===== Set axis limits: =====
xmin_new = -0.5; 
xmax_new = 1.5;
zmin = -0.1; 
zmax = 4;

% Set flag: double check simulations develop AML if their initial starting-position in the state-space lies between c3 and c2 i.e. c3 <= x0 <=c2
maxNumReruns = 5;
attempt = 0;
rerun = true;

rerun_condition = @(x_w0,x0) ((x_w0(1) >= c2_mRNA && x_w0(2) >= c2_mRNA) || (x_w0(1) >= x0(1) && x_w0(2) >= x0(2)));

% Rerun loop:
while rerun && attempt < maxNumReruns
    attempt = attempt + 1;

    % Reset the starting position:
    x = zeros(2, length(t));
    x(:, 1) = x0;

    % Reset storing arrays:
    A11 = []; A12 = []; A21 = [];A22 = [];
    B11 = []; B12 = []; B21 = [];B22 = [];
    Alpha_array = [];
    treatment_time = [];
    effect_matrix = [];

    % ---- Video setup: single combined video ----
    if exist('v', 'var')
        close(v);
    end

    v = VideoWriter('Video1.avi'); %
    v.FrameRate = 10; % MATLAB default is a frame rate of 30 frames per second (fps); Set to a lower frame rate 
    v.Quality   = 95;
    open(v);

    % To fit in a PPT slide
    figW = 768;  % width in pixels = (desired width) x (DPI)
    figH = 384;  % height in pixels = (desired width) x (DPI)
    fig = figure('Units','pixels','Position',[100 100 figW figH],'Color','w','Renderer','opengl');

    % Tiled layout for two panels side-by-side
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax1 = nexttile(tl,1);   % left: 2D trajectory on contour
    ax2 = nexttile(tl,2);   % right: 3D potential

    % Fix limits to avoid flicker across frames
    xlim(ax1,[xmin_new xmax_new]); ylim(ax1,[xmin_new xmax_new]);
    zmin = -0.1; 
    zmax = 2;

    %  ---- Simulation loop  ----
    for i = 2:length(t)
        
        % --- No treatment ---
        if t(i) <= w0 
            Drift = SDE(t(i-1), x(:, i-1));
            F = F_untreated;

         % --- Treatment ---
        else
            % Compute chemo effect at current time point
            effect = rr.*E(find((stime<=t(i)),1,'last'));
            effect_matrix = [effect_matrix;effect];

            norm_effect = effect/effect_max;

            % Update tensor elements
            psi(1,1) = round((psi_untreated(1,1) + (psi_maxEffect(1,1) - psi_untreated(1,1))*norm_effect),2);
            psi(1,2)= round((psi_untreated(1,2) + (psi_maxEffect(1,2) - psi_untreated(1,2))*norm_effect),2);
            psi(2,1) = round((psi_untreated(2,1) + (psi_maxEffect(2,1) - psi_untreated(2,1))*norm_effect),2);
            psi(2,2)= round((psi_untreated(2,2) + (psi_maxEffect(2,2) - psi_untreated(2,2))*norm_effect),2);

            phi(1,1) = round((phi_untreated(1,1) + (phi_maxEffect(1,1) - phi_untreated(1,1))*norm_effect),2);
            phi(1,2)= round((phi_untreated(1,2) + (phi_maxEffect(1,2) - phi_untreated(1,2))*norm_effect),2);
            phi(2,1) = round((phi_untreated(2,1) + (phi_maxEffect(2,1) - phi_untreated(2,1))*norm_effect),2);
            phi(2,2)= round((phi_untreated(2,2) + (phi_maxEffect(2,2) - phi_untreated(2,2))*norm_effect),2);

            % During max effect alpha needs to be alpha_treated and during untreated conditions, alpha = 1
            alpha = round((alpha_treated + (alpha_untreated - alpha_treated)*(alpha_untreated - norm_effect)),2);
            
            % Store information
            A11 = [A11; psi(1,1)]; A12 = [A12;psi(1,2)]; A21 = [A21;psi(2,1)];A22 = [A22;psi(2,2)];
            B11 = [B11;phi(1,1)]; B12 = [B12;phi(1,2)]; B21 = [B21;phi(2,1)];B22 = [B22;phi(2,2)];
            Alpha_array = [Alpha_array;alpha];
            treatment_time = [treatment_time; t(i)];

            % ===== Update drift with corresponding tensor elements =====
            Drift = SDE_chemo(t(i-1), x(:, i-1)); 

            % ===== Update potential =====
            F = alpha.*(psi(1,1).*X.^2 + (psi(1,2) + psi(2,1)).*X.*Y + psi(2,2).*Y.^2 + effect).*((phi(1,1).*(X-c1).^2 + (phi(1,2) + phi(2,1)).*(X-c1).*(Y-c1) + phi(2,2).*(Y-c1).^2) + gamma);

            % ===== store info to plot F at maximum time of chemo effect: =====
            if(isequal(psi,psi_maxEffect) && isequal(phi,phi_maxEffect))
                F_max_effect = F;
                max_effect = effect;
            end
        end
        
        % ===== Brownian noise =====
        dB = sqrt(dt_sim)*randn(2,1);

        % ===== Update particle position according to Langevin equation
        x(:, i) = x(:, i-1) + Drift*dt_sim + Diff*dB;

        % ---------- LEFT SUBPLOT: 2D trajectory ----------
        cla(ax1); hold(ax1,'on');
        contour(ax1, X, Y, log10(F), 30, 'LineWidth',1);
        xticks(ax1, [c3 c2_mRNA c1]); xticklabels(ax1, {'$c_3$','$c_2$','$c_1$'});
        yticks(ax1, [c3 c2_miRNA c1]); yticklabels(ax1, {'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'});
        ax1.XAxis.TickLabelInterpreter = 'latex';
        ax1.YAxis.TickLabelInterpreter = 'latex';
        %Plot trajectory up to current time point
        scatter(ax1,x(1,1),x(2,1),'k','filled'); % start
        plot(ax1,x(1,1:i),x(2,1:i),'k','LineWidth',1);
        if i == length(t)
            scatter(ax1,x(1,i),x(2,i),'c','filled'); 
        end
        set(ax1, 'FontSize', 12);  % sets tick labels for left panel
        xlabel(ax1,'mRNA','FontSize',14);
        ylabel(ax1,'miRNA','FontSize',14);
        title(ax1,'2D simulation','FontSize',16);
        colorbar(ax1,'off'); caxis(ax1,'auto');

        % ---------- RIGHT SUBPLOT: 3D potential ----------
        cla(ax2);
        h = surfc(ax2,X,Y,F);
        hold(ax2,'on');    
        % Set z-limits for axis
        set(ax2,'ZLim',[zmin,zmax]);
        colorbar(ax2,'off'); 
        caxis(ax2,[zmin,zmax]);  
        % Contour levels under the surface
        hContour = h(2);
        nLevels = 30; 
        hContour.LevelList = linspace(zmin,zmax,nLevels);
        % Label x and y axis with critical points
        xticks(ax2,[c3 c2_mRNA c1]); xticklabels(ax2,{'$c_3$','$c_2$','$c_1$'});
        yticks(ax2,[c3 c2_miRNA c1]); yticklabels(ax2,{'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'});
        ax2.XAxis.TickLabelInterpreter = 'latex';
        ax2.YAxis.TickLabelInterpreter = 'latex';
        set(ax2, 'FontSize', 12);  % sets tick labels for left panel
        xlabel(ax2,'mRNA','FontSize',14);
        ylabel(ax2,'miRNA','FontSize',14);
        title(ax2,'Potential','FontSize',16);
        view(ax2,45,18);

        % ---------- Capture combined frame & write ----------
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);

    end
    close(v);

    % --------------Check condition at t= w0---------
    idx_w0 = find(t==3);
    x_w0 = x(:,idx_w0);

    if rerun_condition(x_w0,x0)
        fprintf('DID NOT DEVELOP AML!! RERUNING....');
        close(fig);
        rerun = true;
    else
        rerun = false; %AML developed by t = w0, keep the simulation
    end
end

if attempt == maxNumReruns
    warning('Maximum reruns reached- last simulation kept!')
end

%% ------------------Figures--------------------------

%% Double-well potential: NO treatment (Figure 2B, top panel) 

fig = figure;
ax = gca;
h = surfc(X,Y,F_untreated);
hold on
% Set z-limits for axis
zmin = -0.1;
zmax = 1;
set(ax,'zlim',[zmin,zmax]);
zticks(0:0.2:zmax);
caxis([zmin, zmax]);
% Add contour lines based on z-lim:
hContour = h(2);
nLevels = 30; %number of contour lines
hContour.LevelList = linspace(zmin, zmax, nLevels);
% Label x and y axis with critical points
xticks([c3 c2_mRNA c1]); 
xticklabels({'$c_3$','$c_2$','$c_1$'});
yticks([c3 c2_miRNA c1]);
yticklabels({'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'});
set(gca,'TickLabelInterpreter','latex')
xlabel('mRNA','FontName', 'Arial', 'FontSize', 8, 'FontUnits', 'points');
ylabel('miRNA','FontName', 'Arial', 'FontSize', 8, 'FontUnits', 'points');
view(45,18); 
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'FontUnits', 'points');
set(fig, 'Resize', 'off');
print(fig,'untreated_potential_Fig2B.svg','-dsvg')
exportgraphics(fig,'untreated_potential_Fig2B.png','Resolution',600,'BackgroundColor','none')

%% Untreated potential overlaid with data (Figure 2B, bottom panel) 
fig = figure;
contour(X,Y,log10(F_untreated),30);
hold on;
%data:
PC2_2018CM_mRNA = data.AML2018.mRNA.group_pc_val;
PC1_2018CM_miRNA = data.AML2018.miRNA.group_pc_val;
scatter(PC2_2018CM_mRNA,PC1_2018CM_miRNA,'k',"filled",'DisplayName','Untreated'); hold on;
xlabel('mRNA');
ylabel('miRNA');
xlim([-0.5 1.5]);
ylim([-0.5 1.5]);
set(gca,'xtick',[c3 c2_mRNA c1],'XTickLabel',{'$c_3$','$c_2$','$c_1$'},'TickLabelInterpreter','latex')
set(gca,'ytick',[c3 c2_miRNA c1],'YTickLabel',{'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'},'TickLabelInterpreter','latex')
box on;
hold off;
caxis
print(fig,'untreatedContour_overlaid_data_Fig2B.svg','-dsvg')
exportgraphics(fig,'untreatedContour_overlaid_data_Fig2B.png','Resolution',600,'BackgroundColor','none')

%% Double-well potential: maximum treatment effect (Figure 2C, top panel) 
fig = figure;
ax = gca;
h = surfc(X,Y,F_max_effect);
hold on
% Set z-limits for axis
zmin = -0.1;
zmax = 1;
set(gca,'zlim',[zmin,zmax]);
zticks(0:0.2:zmax);
caxis([zmin, zmax]);
% Add contour lines based on z-lim:
hContour = h(2);
nLevels = 30; %number of contour lines
hContour.LevelList = linspace(zmin, zmax, nLevels);
% Label x and y axis with critical points
xticks([c3 c2_mRNA c1]); 
xticklabels({'$c_3$','$c_2$','$c_1$'});
yticks([c3 c2_miRNA c1]);
yticklabels({'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'});
set(gca,'TickLabelInterpreter','latex')
xlabel('mRNA','FontName', 'Arial', 'FontSize', 8, 'FontUnits', 'points');
ylabel('miRNA','FontName', 'Arial', 'FontSize', 8, 'FontUnits', 'points');
view(45,18); 
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'FontUnits', 'points');
set(fig, 'Resize', 'off');
print(fig,'treated_potential_Fig2C.svg','-dsvg')
exportgraphics(fig,'treated_potential_Fig2C.png','Resolution',600,'BackgroundColor','none')

%% Overlay treated potential with data (Figure 2C, bottom panel) 
fig = figure;
contour(X,Y,log10(F_max_effect),30); hold on;
mRNA_postchemo_pc = data.mRNA.postchemo_pc_val;
miRNA_postchemo_pc = data.miRNA.postchemo_pc_val;
scatter(mRNA_postchemo_pc,miRNA_postchemo_pc,'m',"filled"); hold on;
xticks([c3 c2_mRNA c1]); 
xticklabels({'$c_3$','$c_2$','$c_1$'});
yticks([c3 c2_miRNA c1]);
yticklabels({'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'});
set(gca,'TickLabelInterpreter','latex')
xlim([xmin_new xmax_new])
ylim([xmin_new xmax_new])
hold on
xlabel('mRNA');
ylabel('miRNA')
caxis
print(fig,'treated_contours_data_Fig2C.svg','-dsvg')
exportgraphics(fig,'treated_contours_data_Fig2C.png','Resolution',600,'BackgroundColor','none')

%% mRNA-miRNA trajectory & deformation of potential under maximum treatment effect (Figure 2D) 
fig = figure;
contour(X,Y,log10(F_max_effect),30)
xticks([c3 c2_mRNA c1]); 
xticklabels({'$c_3$','$c_2$','$c_1$'});
yticks([c3 c2_miRNA c1]);
yticklabels({'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'});
set(gca,'TickLabelInterpreter','latex')
xlim([xmin_new xmax_new])
ylim([xmin_new xmax_new])
hold on
plot(x(1,:),x(2,:), 'Color', [137 64 151]/255);%mRNA VS miRNA
scatter(x(1,1),x(2,1), 'k', 'filled'); %start
scatter(x(1,end),x(2,end), 'c','filled'); %end
xlabel('mRNA');
ylabel('miRNA');
title('2D state-transition model');
caxis
print(fig,'2D_model_trajectory_Fig2D.svg','-dsvg')
exportgraphics(fig,'2D_model_trajectory_Fig2D.png','Resolution',600,'BackgroundColor','none')

%% Plot Individual Trajectories (Figure 2E) 

t_shifted = t - w0; %shift time to start at -3

fig = figure;
plot(t_shifted, x(1, :),'Color', [0.3059, 0.6549, 0.1804],'LineWidth',1.5);
hold on
plot(t_shifted, x(2, :),'Color',[0.1294, 0.3725, 0.6039],'LineWidth',1.5);

fill([c3 c1 c1 c3], [xmin_new xmin_new xmax_new xmax_new], ...
     'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlim([-3 t_shifted(end)])
xlabel('Time (weeks)');
ylabel('State-space');
legend('mRNA', 'miRNA');
legend boxoff;
set(gca,'ylim',[xmin_new xmax_new],'YTick',[c3 c2_mRNA c1],'YTickLabel',{'$c_3$','$c_2$','$c_1$'},'TickLabelInterpreter','latex')
title('State-transition trajectories');
print(fig,'single_trajectories_Fig2E.svg','-dsvg')
exportgraphics(fig,'single_trajectories_Fig2E.png','Resolution',600,'BackgroundColor','none')

%% Identify the minima, maxima and/or saddle points
% --- Symbolic setup ---
syms x y real
F = (x.^2 + y.^2).*(((x-c1).^2 + (y-c1).^2) + gamma);
finding_MinMaxSaddle_2D;