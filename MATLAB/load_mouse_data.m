function data = load_mouse_data(params)

    % ================== LOAD DATA ==================
    amldata = readtable('m_miRNA_combined_all.csv');
    amldata.miRNA_PC1 = amldata.miRNA_PC1*(-1);% -amldata.miRNA_PC1;
    
    % ================== COMMON SETTINGS ==================
    time_window = [-3, 10]; % corresponds to the min (t=-3) and max (t=10) experimental times
    max_time_pts = 12;
       
    % ================== RESCALE CRITICAL POINTS & STATE-SPACE TO BE BETWEEN 0-1 ==================
    
    % ---- Rescale mRNA critical points to 0-1 ----
    mRNA_rescale = scale_func(params.mRNA.xmin, params.mRNA.xmax, params.mRNA.raw_cp); %function to scale the state-space
    c_mRNA = apply_scaling(params.mRNA.raw_cp, params.mRNA.xmin, params.mRNA.xmax);
    
    % ---- Rescale miRNA critical points to 0-1----
    miRNA_rescale = scale_func(params.miRNA.xmin, params.miRNA.xmax, params.miRNA.raw_cp); %function to scale the state-space
    c_miRNA = apply_scaling(params.miRNA.raw_cp, params.miRNA.xmin, params.miRNA.xmax);
    
    % ================== MOUSE IDS ==================
    mid_all = [4443;4436;4433;4419;4329;4324;4321;4535;4506];
    rmIDs = [4321 4324]; %remove mice because cKit% < 20%
    keepIdx = ~ismember(mid_all, rmIDs);
    mid = mid_all(~ismember(mid_all, rmIDs));  
    
    % ================== CKIT ==================
    data.ckit = extract_data(amldata, mid,"mouse_id_x", "sample_weeks_x", "percent_ckit_x", @(x)x, time_window, max_time_pts);
    data.ckit.ids = mid;
    data.ckit.keepIdx = keepIdx;

    % ================== mRNA ==================
    data.mRNA = extract_chemo_data(amldata, mid, "mouse_id_x", "sample_weeks_x", "mRNA_PC2", mRNA_rescale, time_window, max_time_pts);
    data.mRNA.c = c_mRNA;
    data.mRNA.ids = mid;
    data.mRNA.keepIdx = keepIdx;
    
    % ================== miRNA ==================
    data.miRNA = extract_chemo_data(amldata, mid,"mouse_id_y", "sample_weeks_y", "miRNA_PC1", miRNA_rescale, time_window, max_time_pts);
    data.miRNA.c = c_miRNA;
    data.miRNA.ids = mid;
    data.miRNA.keepIdx = keepIdx;

    % ================== 2018 AML ==================
    data.AML2018.miRNA = extract_group_data(amldata, "18_CM","mouse_id_y", "sample_weeks_y", "miRNA_PC1", miRNA_rescale);

    data.AML2018.mRNA = extract_group_data(amldata, "18_CM", "mouse_id_x", "sample_weeks_x", "mRNA_PC2", mRNA_rescale);

    % ================== 2018 CONTROL ==================
    data.ctrl2018.miRNA = extract_group_data(amldata, "18_Ctrl","mouse_id_y", "sample_weeks_y", "miRNA_PC1",miRNA_rescale);

    data.ctrl2018.mRNA = extract_group_data(amldata, "18_Ctrl","mouse_id_x", "sample_weeks_x", "mRNA_PC2", mRNA_rescale);
end

