%% Simulation parameters

% ===== Time parameters =====
dt_sim = 0.1;
t = 0:dt_sim:T;
w0 = 3; %start of treatment

% ===== Diffusion coefficient from MSD Analysis (see compute_msd.m): =====
Dx_mRNA = 0.0287;
Dx_miRNA = 0.0244;
Diff = diag([sqrt(2*Dx_mRNA), sqrt(2*Dx_miRNA)]);

% ===== Defining the meshgrid =====
xmin = -1; 
xmax = 2;
myx = linspace(xmin,xmax);
myy = myx;
[X,Y] = meshgrid(myx,myy);

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
chemo_effect_different_doses; 
rr = 1;   %turn chemo effect on/off

%% Define untreated and treated potentials where SDE = [-dF/dx; -dF/dy]

% ===== Untreated potential & drift =====
% Untreated:F(X,Y) = (X^2 + Y^2)*(((X-c1)^2 + (Y-c1)^2) + gamma)
SDE = @(t, x) [-((2.*x(1).*(((x(1)-c1).^2 + (x(2)-c1).^2) + gamma)) + ((x(1).^2 + x(2).^2) .* (2.*(x(1) - c1)))); -((2.*x(2).*(((x(1)-c1).^2 + (x(2)-c1).^2)+ gamma)) + (2.*(x(2)-c1).*(x(1).^2 + x(2).^2)))];

% ===== Treated potential & drift =====
SDE_chemo = @(t,x,effect) [-alpha.*(((2*psi(1,1).*x(1) + (psi(1,2) + psi(2,1)).*x(2)).*((phi(1,1).*(x(1)-c1).^2 + (phi(1,2) + phi(2,1)).*(x(1)-c1).*(x(2)-c1) + phi(2,2).*(x(2)-c1).^2) + gamma)) + (psi(1,1).*x(1).^2 + (psi(1,2) + psi(2,1)).*x(1).*x(2) + psi(2,2).*x(2).^2 + rr.*effect).*(2*phi(1,1).*(x(1)-c1) + (phi(1,2) + phi(2,1)).*(x(2)-c1))); -alpha.*((((psi(1,2) + psi(2,1)).*x(1) + 2*psi(2,2).*x(2)).*((phi(1,1).*(x(1)-c1).^2 + (phi(1,2) + phi(2,1)).*(x(1)-c1).*(x(2)-c1) + phi(2,2).*(x(2)-c1).^2) + gamma)) + ((psi(1,1).*x(1).^2 + (psi(1,2) + psi(2,1)).*x(1).*x(2) + psi(2,2).*x(2).^2 + rr.*effect).*((phi(1,2) + phi(2,1)).*(x(1)-c1) + 2*phi(2,2).*(x(2)-c1))))];

%% Simulate SDE
x0_statespace = zeros(2, length(x0_mRNA));
x0_statespace(1,:) = x0_mRNA;
x0_statespace(2,:) = x0_miRNA;

x = zeros(2, length(t));

%store the final position and if it reached a curative state
x_final = zeros(2,length(dose_factors));
cured_mRNA = zeros(length(x0_statespace),length(dose_factors));
cured_miRNA = zeros(length(x0_statespace),length(dose_factors));
cured_both = zeros(length(x0_statespace),length(dose_factors));

for j = 1:length(dose_factors)

    % Get the Effect (E) associated with the current dose factor gamma_d
    gamma_d = dose_factors(j);
    E = E_matrix(j,:);
    [effect_max, idx_max] = max(E); %maximum effect

    for k = 1:length(x0_statespace)

        % Get the current position
        x0 = x0_statespace(:,k);
        x(:,1) = x0;
       
        % Update the state-space position
        for i = 2:length(t)
            %Chemo:
            if t(i) <= w0
              Drift = SDE(t(i-1), x(:, i-1));
              
              F = (X.^2 + Y.^2).*(((X-c1).^2 + (Y-c1).^2) + gamma);

            else
                % Compute chemo effect at current time
                effect = rr.*E(find((stime<=t(i)),1,'last')); %Effect at the current time point

                norm_effect = effect/effect_max;

                % tensor elements
                psi(1,1) = round((psi_untreated(1,1) + (psi_maxEffect(1,1) - psi_untreated(1,1))*norm_effect),2);
                psi(1,2) = round((psi_untreated(1,2) + (psi_maxEffect(1,2) - psi_untreated(1,2))*norm_effect),2);
                psi(2,1) = round((psi_untreated(2,1) + (psi_maxEffect(2,1) - psi_untreated(2,1))*norm_effect),2);
                psi(2,2) = round((psi_untreated(2,2) + (psi_maxEffect(2,2) - psi_untreated(2,2))*norm_effect),2);

                phi(1,1) = round((phi_untreated(1,1) + (phi_maxEffect(1,1) - phi_untreated(1,1))*norm_effect),2);
                phi(1,2) = round((phi_untreated(1,2) + (phi_maxEffect(1,2) - phi_untreated(1,2))*norm_effect),2);
                phi(2,1) = round((phi_untreated(2,1) + (phi_maxEffect(2,1) - phi_untreated(2,1))*norm_effect),2);
                phi(2,2) = round((phi_untreated(2,2) + (phi_maxEffect(2,2) - phi_untreated(2,2))*norm_effect),2);

                %During max effect alpha needs to be alpha_treated and during untreated conditions, alpha = 1
                alpha = round((alpha_treated + (alpha_untreated - alpha_treated)*(alpha_untreated - norm_effect)),2);

                % Update drift with corresponding tensor elements
                Drift = SDE_chemo(t(i-1), x(:, i-1),effect); %apply chemo

                % Update F
                F = alpha.*(psi(1,1).*X.^2 + (psi(1,2) + psi(2,1)).*X.*Y + psi(2,2).*Y.^2 + effect).*((phi(1,1).*(X-c1).^2 + (phi(1,2) + phi(2,1)).*(X-c1).*(Y-c1) + phi(2,2).*(Y-c1).^2) + gamma);

                %store info to plot F at maximum time of chemo effect:
                if(psi(1,1)==psi_maxEffect(1,1) && psi(1,2)==psi_maxEffect(1,2) && psi(2,1)==psi_maxEffect(2,1) && psi(2,2)==psi_maxEffect(2,2) && phi(1,1)==phi_maxEffect(1,1) && phi(1,2)==phi_maxEffect(1,2) && phi(2,1)==phi_maxEffect(2,1) && phi(2,2)==phi_maxEffect(2,2))
                    F_max_effect = F;
                    max_effect = effect;
                end
            end

            dB = sqrt(dt_sim)*randn(2,1); %Brownian inc

            x(:, i) = x(:, i-1) + Drift*dt_sim + Diff*dB; %update position            
        end

        % ADD A FLAG: some trajectories do not reach AML state and immediately jump to a health state (innacurate). If this occurs need to rerun sim
        idx_w0 = find(t==3);
        x_w0 = x(:,idx_w0);

        if ((x0(1,:)<c2_mRNA && x0(2,:)<c2_mRNA) && (x0(1,:)>c3 && x0(2,:)>c3))
            %Apply flag only for starting position below c2
            
            if((x_w0(1,:) >= c2_mRNA && x_w0(2,:) >= c2_mRNA) || (x_w0(1,:) >= x0(1,:) && x_w0(2,:) >= x0(2,:)))
                disp('DID NOT DEVELOP AML!! RERUNING....')
                disp(gamma_d);
                disp(x0);
    
                % Rerun:
                for i = 2:length(t)
                    %Chemo:
                    if t(i) <= w0
                      Drift = SDE(t(i-1), x(:, i-1));
                     
                      F = (X.^2 + Y.^2).*(((X-c1).^2 + (Y-c1).^2) + gamma);
                    else
                        % Compute chemo effect at current time
                        effect = rr.*E(find((stime<=t(i)),1,'last')); %Effect at the current time point
    
                        norm_effect = effect/effect_max;
    
                        % tensor elements
                        psi(1,1) = round((psi_untreated(1,1) + (psi_maxEffect(1,1) - psi_untreated(1,1))*norm_effect),2);
                        psi(1,2) = round((psi_untreated(1,2) + (psi_maxEffect(1,2) - psi_untreated(1,2))*norm_effect),2);
                        psi(2,1) = round((psi_untreated(2,1) + (psi_maxEffect(2,1) - psi_untreated(2,1))*norm_effect),2);
                        psi(2,2) = round((psi_untreated(2,2) + (psi_maxEffect(2,2) - psi_untreated(2,2))*norm_effect),2);
    
                        phi(1,1) = round((phi_untreated(1,1) + (phi_maxEffect(1,1) - phi_untreated(1,1))*norm_effect),2);
                        phi(1,2) = round((phi_untreated(1,2) + (phi_maxEffect(1,2) - phi_untreated(1,2))*norm_effect),2);
                        phi(2,1) = round((phi_untreated(2,1) + (phi_maxEffect(2,1) - phi_untreated(2,1))*norm_effect),2);
                        phi(2,2) = round((phi_untreated(2,2) + (phi_maxEffect(2,2) - phi_untreated(2,2))*norm_effect),2);
    
                        %During max effect alpha needs to be alpha_treated and during untreated conditions, alpha = 1
                        alpha = round((alpha_treated + (alpha_untreated - alpha_treated)*(alpha_untreated - norm_effect)),2);
    
                        % Update drift with corresponding tensor elements
                        Drift = SDE_chemo(t(i-1), x(:, i-1),effect); %apply chemo
                           
                        % Update F
                        F = alpha.*(psi(1,1).*X.^2 + (psi(1,2) + psi(2,1)).*X.*Y + psi(2,2).*Y.^2 + effect).*((phi(1,1).*(X-c1).^2 + (phi(1,2) + phi(2,1)).*(X-c1).*(Y-c1) + phi(2,2).*(Y-c1).^2) + gamma);
    
                        %store info to plot F at maximum time of chemo effect:
                        if(psi(1,1)==psi_maxEffect(1,1) && psi(1,2)==psi_maxEffect(1,2) && psi(2,1)==psi_maxEffect(2,1) && psi(2,2)==psi_maxEffect(2,2) && phi(1,1)==phi_maxEffect(1,1) && phi(1,2)==phi_maxEffect(1,2) && phi(2,1)==phi_maxEffect(2,1) && phi(2,2)==phi_maxEffect(2,2))
                            F_max_effect = F;
                            max_effect = effect;
                        end
                    end
    
                    dB = sqrt(dt_sim)*randn(2,1); %Brownian inc
    
                    x(:, i) = x(:, i-1) + Drift*dt_sim + Diff*dB; %update position
    
                end
    
                %Double check the issue is fixed
                idx_w0 = find(t==3);
                x_w0 = x(:,idx_w0);
                if((x_w0(1,:) >= c2_mRNA && x_w0(2,:) >= c2_mRNA) || (x_w0(1,:) >= x0(1,:) && x_w0(2,:) >= x0(2,:)))
                    disp('NOT FIXED!!!!!! Rerunning again...')

                    % Rerun:
                    for i = 2:length(t)
                        %Chemo:
                        if t(i) <= w0
                          Drift = SDE(t(i-1), x(:, i-1));
                         
                          F = (X.^2 + Y.^2).*(((X-c1).^2 + (Y-c1).^2) + gamma);
                        else
                            % Compute chemo effect at current time
                            effect = rr.*E(find((stime<=t(i)),1,'last')); %Effect at the current time point
        
                            norm_effect = effect/effect_max;
        
                            % tensor elements
                            psi(1,1) = round((psi_untreated(1,1) + (psi_maxEffect(1,1) - psi_untreated(1,1))*norm_effect),2);
                            psi(1,2) = round((psi_untreated(1,2) + (psi_maxEffect(1,2) - psi_untreated(1,2))*norm_effect),2);
                            psi(2,1) = round((psi_untreated(2,1) + (psi_maxEffect(2,1) - psi_untreated(2,1))*norm_effect),2);
                            psi(2,2) = round((psi_untreated(2,2) + (psi_maxEffect(2,2) - psi_untreated(2,2))*norm_effect),2);
        
                            phi(1,1) = round((phi_untreated(1,1) + (phi_maxEffect(1,1) - phi_untreated(1,1))*norm_effect),2);
                            phi(1,2) = round((phi_untreated(1,2) + (phi_maxEffect(1,2) - phi_untreated(1,2))*norm_effect),2);
                            phi(2,1) = round((phi_untreated(2,1) + (phi_maxEffect(2,1) - phi_untreated(2,1))*norm_effect),2);
                            phi(2,2) = round((phi_untreated(2,2) + (phi_maxEffect(2,2) - phi_untreated(2,2))*norm_effect),2);
        
                            %During max effect alpha needs to be alpha_treated and during untreated conditions, alpha = 1
                            alpha = round((alpha_treated + (alpha_untreated - alpha_treated)*(alpha_untreated - norm_effect)),2);
        
                            % Update drift with corresponding tensor elements
                            Drift = SDE_chemo(t(i-1), x(:, i-1),effect); %apply chemo
        
                            % Update F
                            F = alpha.*(psi(1,1).*X.^2 + (psi(1,2) + psi(2,1)).*X.*Y + psi(2,2).*Y.^2 + effect).*((phi(1,1).*(X-c1).^2 + (phi(1,2) + phi(2,1)).*(X-c1).*(Y-c1) + phi(2,2).*(Y-c1).^2) + gamma);

                            %store info to plot F at maximum time of chemo effect:
                            if(psi(1,1)==psi_maxEffect(1,1) && psi(1,2)==psi_maxEffect(1,2) && psi(2,1)==psi_maxEffect(2,1) && psi(2,2)==psi_maxEffect(2,2) && phi(1,1)==phi_maxEffect(1,1) && phi(1,2)==phi_maxEffect(1,2) && phi(2,1)==phi_maxEffect(2,1) && phi(2,2)==phi_maxEffect(2,2))
                                F_max_effect = F;
                                max_effect = effect;
                            end
                        end
        
                        dB = sqrt(dt_sim)*randn(2,1); %Brownian inc
        
                        x(:, i) = x(:, i-1) + Drift*dt_sim + Diff*dB; %update position
                    end

                    %Double check the issue is fixed
                    idx_w0 = find(t==3);
                    x_w0 = x(:,idx_w0);
                    if((x_w0(1,:) >= c2_mRNA && x_w0(2,:) >= c2_mRNA) || (x_w0(1,:) >= x0(1,:) && x_w0(2,:) >= x0(2,:)))
                        disp('NOT FIXED!!!!!!')
                    else
                        disp('Fixed!')
                    end

                else
                    disp('Fixed!')
                end
    
            end
        end
        
        %Store the final position of each trajectory
        x_final(:,j) = x(:,end);

        % If the last position is greater than c2, the trajectories reach a curative state
        if x(1,end)>c2_mRNA
            cured_mRNA(k,j) = 1;
        else
            cured_mRNA(k,j) = 0;
        end

        if x(2,end)>c2_miRNA
            cured_miRNA(k,j) = 1;
        else
            cured_miRNA(k,j) = 0; 
        end

        if (x(1,end)>c2_mRNA && x(2,end)>c2_miRNA)
            cured_both(k,j) = 1;
        else
            cured_both(k,j) = 0; 
        end
    end
end
