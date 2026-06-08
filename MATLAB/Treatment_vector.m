%% ==============================================================================================
% Largest change in direction between F and F_Rx: Worst directional mismatch between ∇F and ∇F_Rx
% c_min = min((∇F*∇F_Rx)/ (||∇F|| ||∇F_Rx||))
% theta_max = arccos(c_min)
% ===============================================================================================

%% Defining the meshgrid/domain

xmin = -1; 
xmax = 2; 
myx = linspace(xmin,xmax);
myy = myx;
[X,Y] = meshgrid(myx,myy);

%% Parameters to define the potential energy landscape at max treatment effect:
effect_max = max(E);
alpha = alpha_treated;

%Tensor elements at max treatment effect:
psi11 = psi_maxEffect(1,1);
psi12 = psi_maxEffect(1,2);
psi21 = psi_maxEffect(2,1);
psi22 = psi_maxEffect(2,2);

phi11 = phi_maxEffect(1,1);
phi12 = phi_maxEffect(1,2);
phi21 = phi_maxEffect(2,1);
phi22 = phi_maxEffect(2,2);

%% Define the potentials
 
F = F_untreated; 
F_chemo = alpha.*(psi11.*X.^2 + (psi12 + psi21).*X.*Y + psi22.*Y.^2 + effect_max).*((phi11.*(X-c1).^2 + (phi12 + phi21).*(X-c1).*(Y-c1) + phi22.*(Y-c1).^2) + gamma);

%% Compute the force vectors i.e. the negative gradients

% --------------------- No Treatment ---------------------
[dFX, dFY] = gradient(F);
dFX = -dFX;
dFY = -dFY;

% --------------------- Chemotherapy ---------------------
[dFX_chemo, dFY_chemo] = gradient(F_chemo);
dFX_chemo = -dFX_chemo;
dFY_chemo = -dFY_chemo;

%% Compute the norm and dot product 

normdF = sqrt(dFX.^2 + dFY.^2); %no treatment 

normdF_Rx = sqrt(dFX_chemo.^2 + dFY_chemo.^2); %chemotherapy

dotProd = dFX .* dFX_chemo + dFY .* dFY_chemo; %dot product

%% Cosine of angle between gradients
cosTheta = dotProd ./ (normdF .* normdF_Rx);

% remove undefined points (zero gradients)
valid = (normdF > 1e-10) & (normdF_Rx  > 1e-10);
cosTheta(~valid) = NaN;

% numerical safety to keep cosTheta in [-1,1]
cosTheta = max(min(cosTheta,1),-1);

%% Compute c_min (value & location) - biggest direction disagreement
[c_min_all, min_idx_all] = min(cosTheta,[], 'omitnan'); %row vector containing min value of each col of cosTheta
[c_min, min_idx] = min(cosTheta(:),[], 'omitnan');
[row, col] = ind2sub(size(X), min_idx);

% spatial location:
x_min = X(row,col);
y_min = Y(row,col);

% maximum angular difference between the two vector fields
theta_max = acos(c_min);

fprintf('c_min = %.4f\n', c_min);
fprintf('theta_max = %.4f rad (%.2f deg)\n', theta_max, rad2deg(theta_max));
fprintf('Largest directional disagreement at (x*,y*) = (%.4f, %.4f)\n', x_min, y_min);

%% A different way to get the maximum angle (step to double check)
theta = acos(cosTheta);
[maxTheta,max_idx] = max(theta(:));
[row_max, col_max] = ind2sub(size(theta), max_idx);

x_max = X(row,col);
y_max = Y(row,col);

fprintf('Maximum angle = %.2f degrees\n',rad2deg(maxTheta));
fprintf('Occurs at:\n');
fprintf('x = %.4f\n',x_max);
fprintf('y = %.4f\n',y_max);

%% =========== Figure: Vizualizing the angle heatmap and max point (Figure S4B) =========== 

theta_values = acos(cosTheta);
theta_values(~valid) = NaN;
theta_values_degrees = rad2deg(theta_values);

xmin_new = -0.5; 
xmax_new = 1.5;

fig = figure;
ax = gca;
set(ax,'FontName', 'Arial', 'FontSize', 8, 'FontUnits', 'points')
imagesc(myx,myy,theta_values_degrees);
axis image;
set(gca, 'YDir','normal');
%axis xy equal tight;
colormap(turbo);
cb = colorbar;
set(cb,'FontName', 'Arial', 'FontSize', 10)
hold on;
plot(x_min,y_min, 'wo', 'MarkerSize', 8, 'LineWidth', 2);
offset = 0.05;
text(x_min+offset,y_min+offset, sprintf(' \\theta_{max} = %.1f^{\\circ}', rad2deg(theta_max)), 'Color', 'w', 'FontWeight','bold', 'FontName','Arial', 'FontSize',10)
title('Angle between -∇U_p and -∇U_{p}^{Rx}');
xticks([c3 c2_mRNA c1]); 
xticklabels({'$c_3$','$c_2$','$c_1$'});
yticks([c3 c2_miRNA c1]);
yticklabels({'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'});
set(gca,'TickLabelInterpreter','latex')
xlim([xmin_new xmax_new])
ylim([xmin_new xmax_new])
xlabel('mRNA','FontName', 'Arial', 'FontSize', 10, 'FontUnits', 'points');
ylabel('miRNA','FontName', 'Arial', 'FontSize', 10, 'FontUnits', 'points');
print(fig,'max_angle_FigS2B.svg','-dsvg')
exportgraphics(fig,'max_angle_FigS2B.png','Resolution',600,'BackgroundColor','none')


%% ****************** Compute the treatment vector R ******************
% R^ = -∇F^_Rx - (-∇F^) evaluated at (x*,y*)

%% Unit vectors at location of worst directional mismatch
udFX = dFX ./ normdF;
udFY = dFY ./ normdF;
udFX_chemo = dFX_chemo ./ normdF_Rx;
udFY_chemo = dFY_chemo ./ normdF_Rx;

% difference vector field
scale = 0.1;
uR_x = scale.*(udFX_chemo - udFX);
uR_y = scale.*(udFY_chemo - udFY);

%% Extract the treatment vector at the maxima-angle location
treatment_vec = [uR_x(row,col); uR_y(row,col)]; %uF_Rx - uF;%uF - uF_Rx;% uF_Rx - uF;
treatment_vec_mag = norm(treatment_vec);
fprintf('Rx Vec = -∇F^_{Rx} - (-∇F^) evaluated at location of (x*,y*) = (%.4f, %.4f)\n', treatment_vec(1), treatment_vec(2)')
fprintf('||-∇F_{Rx} - (-∇F)|| = %.4f\n', treatment_vec_mag);

%% ==================================== Figures ==========================================================
%% Overlay vector fields (Figure S4C)
ds = 1;  % downsampling factor

fig = figure; 
axMain = axes;
hold(axMain, 'on')

% Original field (black)
quiver(axMain,X(1:ds:end,1:ds:end), Y(1:ds:end,1:ds:end), udFX(1:ds:end,1:ds:end), udFY(1:ds:end,1:ds:end), 0.8, 'k');
% Rx field (magenta)
quiver(axMain,X(1:ds:end,1:ds:end), Y(1:ds:end,1:ds:end), udFX_chemo(1:ds:end,1:ds:end), udFY_chemo(1:ds:end,1:ds:end), 0.8, 'm');
plot(axMain,x_min,y_min, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
axis(axMain, 'equal', 'tight');
box(axMain, 'on')
xlabel(axMain,'mRNA'); 
ylabel(axMain,'miRNA');
title(axMain,'-∇U_p vs -∇U_{p}^{Rx}');
legend(axMain,'-∇U_p','-∇U_{p}^{Rx}', 'Location','northeastoutside');
xticks(axMain,[c3 c2_mRNA c1])
xticklabels(axMain,{'$c_3$','$c_2$','$c_1$'})
yticks(axMain,[c3 c2_miRNA c1])
yticklabels(axMain,{'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'})
axMain.XAxis.TickLabelInterpreter = 'latex';
axMain.YAxis.TickLabelInterpreter = 'latex';
xlim(axMain,[xmin_new xmax_new])
ylim(axMain,[xmin_new xmax_new])
%inset subfigure:
axInset = axes('Position',[0.75 0.3 0.22 0.38]); %[x y w h]
hold(axInset,'on')
box(axInset, 'on')
set(axInset, 'LineWidth', 1.5)
quiver(axInset,X(1:ds:end,1:ds:end), Y(1:ds:end,1:ds:end), udFX(1:ds:end,1:ds:end), udFY(1:ds:end,1:ds:end), 0.8, 'k');
quiver(axInset,X(1:ds:end,1:ds:end), Y(1:ds:end,1:ds:end), udFX_chemo(1:ds:end,1:ds:end), udFY_chemo(1:ds:end,1:ds:end), 0.8, 'm');
plot(axInset,x_min,y_min, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
axis(axInset, 'equal')
set(axInset, 'XTick', [], 'YTick',[])
xlim(axInset, [c2_mRNA - 0.18 c2_mRNA-0.03])
ylim(axInset,[c3+0.12 c3 + 0.24])
xl = xlim(axInset);
yl = ylim(axInset);
rectangle(axMain, 'Position',[xl(1) yl(1) diff(xl) diff(yl)], 'EdgeColor','k', 'LineStyle', '-','LineWidth', 1);
print(fig,'vector_field_FigS2C.svg','-dsvg')
exportgraphics(fig,'vector_field_FigS2C.png','Resolution',600,'BackgroundColor','none')

%% Treatment vector (red) and difference field (gray) (Figure 2F) 
ds = 4; %down-sampling factor

fig = figure;
quiver(X(1:ds:end,1:ds:end), Y(1:ds:end,1:ds:end), uR_x(1:ds:end,1:ds:end), uR_y(1:ds:end,1:ds:end),0,'Color',[.6 .6 .6]); hold on;
quiver(x_min, y_min, treatment_vec(1), treatment_vec(2),0, 'r','LineWidth', 2, 'MaxHeadSize', 1);
axis equal tight;
xlabel('mRNA'); ylabel('miRNA');
title('Treatment vector');
xticks([c3 c2_mRNA c1]); 
xticklabels({'$c_3$','$c_2$','$c_1$'});
yticks([c3 c2_miRNA c1]);
yticklabels({'$\tilde{c}_3$','$\tilde{c}_2$','$\tilde{c}_1$'});
set(gca,'TickLabelInterpreter','latex')
xlim([c3-0.2 c1+0.2])
ylim([c3-0.2 c1+0.2])
print(fig,'treatment_vector_Fig2F.svg','-dsvg')
exportgraphics(fig,'treatment_vector_Fig2F.png','Resolution',600,'BackgroundColor','none')

%% Eigenvalues and eigenvectors
% Psi
fprintf('-------EIGENSTRUCTURE-------')
psi = [psi11 psi12; psi21 psi22];
[V_psi, D_psi] = eig(psi)

V_psi_Emax = effect_max*V_psi

% Phi
phi = [phi11 phi12; phi21 phi22];
[V_phi, D_phi] = eig(phi)
V_phi_Emax = effect_max*V_phi
