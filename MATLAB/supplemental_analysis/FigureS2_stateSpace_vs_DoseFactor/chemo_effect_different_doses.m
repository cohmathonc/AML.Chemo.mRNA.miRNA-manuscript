%This file was set up so that we agree with time in data from -3 to 10

%Time:
t0 = 0; %initial time point (weeks)
tff = 10; %time at which we truncate the data
nt = 1000; %number of time steps

%Treatment:  
w0 = 3; %treatment start at 6 weeks
ti = [1 2 3 4 5]/7 + w0 ; %treatment start times (days)
tj = [1 3 5]/7 + w0 ; %treatment end times (days)

% Dose
Dc_base = 50*0.0257; %treatment effect parameter
Dd_base = 1.5*0.0257; %treatment effect parameter

% half-life 40.4 hours (0.2404762 week) for cytarabine,31.5 hours (0.18751 week) daunorubicin
% cytarabine (50mg/kg/day; 5 days) and daunorubicin (1.5mg/kg/day, 3 days)
% mouse average weight 25.7 gram = 0.0257 kg
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314217/
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4520927/#:~:text=Median%20half%2Dlife%20was%2031.1,days%20after%20the%20last%20dose
lam1 = log(2)/0.2404762;
lam2 = log(2)/0.1875;

%Evaluate the identified PK-PD model at nt timepoints
stime = linspace(0,T,nt);

dt = stime(2)-stime(1);

ke_miRNA = 0.42;      %rate constant for drug transfer PK -> PD compartment
ke_mRNA = 0.60;      %rate constant for drug transfer PK -> PD compartment
Em = pi();      %maximum effect
EC50 = 0.24*Em; %concentration at which 50% of max effect is observed

%Store info:
Dc_array = zeros(1,length(dose_factors));
Dd_array = zeros(1,length(dose_factors));
C_matrix = zeros(length(dose_factors),nt);
E_matrix = zeros(length(dose_factors),nt); 

for i = 1:length(dose_factors)
    gamma_d = dose_factors(i);

    %scale the dose and store
    Dc = gamma_d * Dc_base;
    Dd = gamma_d * Dd_base;

    Dc_array(i) = Dc;
    Dd_array(i) = Dd;
    
    %PK input
    C =  Treatment(stime,ti,tj,Dc,Dd,lam1,lam2); %chemo treatment
    C_matrix(i,:) = C;

    %Effect - miRNA
    E_miRNA = getE(Em,EC50,getCe(ke_miRNA,stime,dt,C)); %Effect = (E_max*C_e)/(EC_50 + Ce)

    %Effect - mRNA
    E_mRNA = getE(Em,EC50,getCe(ke_mRNA,stime,dt,C)); %Effect = (E_max*C_e)/(EC_50 + Ce)
    
    % Combined effect
    E = (E_mRNA + E_miRNA)/2;
    E_matrix(i,:) = E;
    
end
