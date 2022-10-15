% initialize_ABCDCon_MGS
clear all; close all

USING_HALLE_DATASET = 1; % RHit + HC
RHIT_FILES = 1; % Rhit
ASHS_files = 1;

disp('loading rhit set of items')
analMRIDir = 'E:\spatcon_hdz\data_p\analyzed_mri_mixed_items\';

% ASHS tracings
filename = fullfile(analMRIDir, '\multivariate_sanityCheck','allData_spaceCond_Rhit_QA_runs_rmv.csv');
alldata_table1 = readtable(filename,'Delimiter',',');
% manual tracings
filename = fullfile(analMRIDir, '\multivariate_sanityCheck','allData_manual_ROI_spaceCond_Rhit_QA_runs_rmv.csv');
alldata_table2 = readtable(filename,'Delimiter',',');

alldata_table = [alldata_table1; alldata_table2];
clear alldata_table2 alldata_table1;

% fisher zscore data
A = alldata_table.r;
alldata_table.z = 1/2*(log((1+A)./(1-A)));
roi_str = unique(alldata_table.roi);

subj = unique(alldata_table.subj);
rows = find(ismember(alldata_table.roi,roi_str{4}));
phc_ant_table = alldata_table(rows,:);
for ii_s = 1:length(subj)
    N(ii_s) = sum(ismember(phc_ant_table.subj,subj{ii_s}));
end

% We remove the sameHouse_eitherRoom condition which is included in other
% Prep an alternative table comparing sameHouse/diffHouse
A = find(ismember(alldata_table.condition, {'sameHouse_eitherRoom','sameVideo'} ));
alldata_table_house = alldata_table; % removing sameHouse_eitherRoom because it overlaps with other options
alldata_table_house(A,:) = [];
unique(alldata_table_house.condition)
% Conditions left are all in diff videos - same house or diff house

% Remove items with wrong encoding during post-viewing testing
A = find(ismember(alldata_table_house.condition1,'unknownRoom'));
alldata_table_house(A,:) = [];

row_num1 = find(ismember(alldata_table_house.entryRoom,'entryRoom'));
alldata_table_house.room_cond(row_num1) = {'sameRoom'};
row_num2 = find(ismember(alldata_table_house.entryRoom,'innerRoom'));
alldata_table_house.room_cond(row_num2) = {'sameRoom'};
rest_rows = ~ismember(1:length(alldata_table_house.room_cond),[row_num1',row_num2']);
alldata_table_house.room_cond(rest_rows) = {'diffRoom'};
unique(alldata_table_house.room_cond)

% Generate furthest space conditions
row_num1 = find(ismember(alldata_table_house.room_cond,'sameRoom') & ...
    ismember(alldata_table_house.condition,'differentVideo_sameHouse'));
alldata_table_house.spaceCond(row_num1) = {'sameHouseRoom_diffVid'};
row_num2 = find(ismember(alldata_table_house.room_cond,'diffRoom') & ...
    ismember(alldata_table_house.condition,'differentHouse_eitherRoom'));
alldata_table_house.spaceCond(row_num2) = {'diffHouseRoom_diffVid'};
rest_rows = find(~ismember(1:length(alldata_table_house.spaceCond),[row_num1',row_num2']));
alldata_table_house.spaceCond(rest_rows) = {'diffHouseSameRoom_diffVid'};
unique(alldata_table_house.spaceCond)

% Generate furthest space conditions - option 2
row_num1 = find(ismember(alldata_table_house.room_cond,'diffRoom') & ...
    ismember(alldata_table_house.condition,'differentHouse_eitherRoom'));
alldata_table_house.spaceCond_joint(row_num1) = {'diffHouseRoom_diffVid'};
rest_rows = find(~ismember(1:length(alldata_table_house.spaceCond_joint),row_num1));
alldata_table_house.spaceCond_joint(rest_rows) = {'other'};
unique(alldata_table_house.spaceCond_joint)


% Left hemisphere only, correctly encoded items
A = find(ismember(alldata_table_house.hemi, 'right' )); % remove right-hem entries
alldata_table_house_left = alldata_table_house;
alldata_table_house_left(A,:) = [];
A = find(ismember(alldata_table_house_left.condition1,'unknownRoom')); % remove items with wrong encoding during testing
alldata_table_house_left(A,:) = [];

% Left hemisphere only, correctly encoded items
A = find(ismember(alldata_table_house.hemi, 'left' )); % remove left-hem entries
alldata_table_house_right = alldata_table_house;
alldata_table_house_right(A,:) = [];
A = find(ismember(alldata_table_house_right.condition1,'unknownRoom')); % remove items with wrong encoding during testing
alldata_table_house_right(A,:) = [];

if false % Exploratory data analysis to test the sensitivity to different aspects of space in the data - proximity to entry door
    
    formula = 'z~spaceCond_joint*roi + spaceCond_joint*hemi + roi*hemi + (1|subj)';
    lme1 = fitlme(alldata_table_house,formula);
    formula = 'z~spaceCond_joint*roi*hemi + (1|subj)';
    lme2 = fitlme(alldata_table_house,formula);
    results = compare(lme1,lme2);
    disp('3-way interaction - spaceCond_joint, roi, hem')
    disp(results.pValue)
    % Significant hem/roi interaction
    % significant interaction P = 0.03
    % For spaceCond:
    % significant interaction P = 0.002
    
    row_num1 = find(ismember(alldata_table_house.roi,'CA1_body'));
    row_num2 = find(ismember(alldata_table_house.roi,'CA2_3_DG_body'));
    alldata_table_house_ca1_ca23dg = alldata_table_house([row_num1(:)', row_num2(:)'],:) ;
    formula = 'z~spaceCond_joint*roi + spaceCond_joint*hemi + roi*hemi + (1|subj)';
    lme1 = fitlme(alldata_table_house_ca1_ca23dg,formula);
    formula = 'z~spaceCond_joint*roi*hemi + (1|subj)';
    lme2 = fitlme(alldata_table_house_ca1_ca23dg,formula);
    results = compare(lme1,lme2);
    disp('3-way interaction - spaceCond_joint, roi, hem')
    disp(results.pValue)
    
    
    formula = 'z~spaceCond_joint + roi + (1|subj)';
    lme1 = fitlme(alldata_table_house_ca1_ca23dg,formula);
    formula = 'z~spaceCond_joint*roi + (1|subj)';
    lme2 = fitlme(alldata_table_house_ca1_ca23dg,formula);
    results = compare(lme1,lme2);
    disp('2-way interaction - spaceCond_joint, roi, in left hem')
    disp(results.pValue)
    % significant interaction < 10^-4
    
    formula = 'z~spaceCond_joint + roi + (1|subj)';
    lme1 = fitlme(alldata_table_house_right,formula);
    formula = 'z~spaceCond_joint*roi + (1|subj)';
    lme2 = fitlme(alldata_table_house_right,formula);
    results = compare(lme1,lme2);
    disp('2-way interaction - spaceCond_joint, roi, in right hem')
    disp(results.pValue)
    % n.s interaction
    
    % Mixed effect model for space condition
    pvalue_vec = [];  test_str = [];
    for hem_i = 3
        for roi_i = 1:7
            if hem_i == 1
                alldata_table_house_hem = alldata_table_house_left;
            elseif hem_i == 2
                alldata_table_house_hem = alldata_table_house_right;
            elseif hem_i == 3
                alldata_table_house_hem = alldata_table_house;
            end
            
            row_num = find(ismember(alldata_table_house_hem.roi,roi_str{roi_i}));
            alldata_table_house_hem_roi = alldata_table_house_hem(row_num,:);
            
            roi_test{roi_i} = unique(alldata_table_house_hem_roi.roi);
            condition_val = unique(alldata_table_house_hem_roi.condition);
            spaceCond_val = unique(alldata_table_house_hem_roi.spaceCond_joint);
            room_cond_val = unique(alldata_table_house_hem_roi.room_cond);
            
            test_i = 1;
            formula = 'z~ 1  + (1|subj)';
            lme1 = fitlme(alldata_table_house_hem_roi,formula);
            formula = 'z~ spaceCond_joint + (1|subj)';
            lme2 = fitlme(alldata_table_house_hem_roi,formula);
            results = compare(lme1,lme2);
            test_str{test_i} = 'DV/DH/DR vs all others';
            pvalue_vec(hem_i, roi_i,test_i) = results.pValue(2);
            
            row_num = find(ismember(alldata_table_house_hem_roi.condition,condition_val{2}));
            alldata_table_samehouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            test_i = 2;
            formula = 'z ~ 1  + (1|subj)';
            lme1 = fitlme(alldata_table_samehouse_hem_roi,formula);
            formula = 'z ~ entryRoom + (1|subj)';
            lme2 = fitlme(alldata_table_samehouse_hem_roi,formula);
            results = compare(lme1,lme2);
            test_str{test_i} = 'DV/SH entryRoom';
            pvalue_vec(hem_i, roi_i,test_i) = results.pValue(2);
            
        end
    end
    
end % various mixed effect models

% 'DV/DH/DR vs all others' = pmERC and alERC are sig.

% Hypothesis - higher sensitivity for space will be correlated with better
% judgement, 'memorability'.


% define an RSA grade per subject for the difference between sameH/sameR
% and diff house items.
%% Add behacior to the story
% subject sepcific onset files are extracted using stimuliExtraction.m
onsetFiles_root = 'E:\spatcon_hdz\data_p\onsetFiles\';
subj_list = unique(alldata_table.subj);
for ii_s = 1:length(subj_list)
    subjects(ii_s) = str2num(subj_list{ii_s}(end-2:end));
end

% add behavioral aspect -
clear HC_RC_perc_all
for subj_i = 1:length(subjects)
    subj = subjects(subj_i);
    
    folderName = fullfile(onsetFiles_root, sprintf('s%03d',subj));
    try
        stimuliFile = fullfile(folderName,sprintf('stimuliList_s%03d.mat',subj));
        mm = matfile(stimuliFile);
        behavStats = mm.behavStats;
    catch
        warning('onset file not found')
    end
    
    HC_RC_perc_all(subj_i) = behavStats.space_correct_perc;
end

% add behavioral aspect -
clear rhit_perc_all
for subj_i = 1:length(subjects)
    subj = subjects(subj_i);
    
    folderName = fullfile(onsetFiles_root, sprintf('s%03d',subj));
    try
        stimuliFile = fullfile(folderName,sprintf('stimuliList_s%03d.mat',subj));
        mm = matfile(stimuliFile);
        behavStats = mm.behavStats;
    catch
        warning('onset file not found')
    end
    
    rhit_perc_all(subj_i) = behavStats.rhit_perc;
end



% Adding RSA spatial score
%%
pvalue_vec = [];  test_str = []; pvalue_ranksum_vec = []; pval_all_subj = [];
clear R_space R_rhit corr_p_space corr_p_rhit RSA_score_SUBJ_ROI
for hem_i = 1:2
    for roi_i = 1:7
        if hem_i == 1
            alldata_table_house_hem = alldata_table_house_left;
        else
            alldata_table_house_hem = alldata_table_house_right;
        end
        
        row_num = find(ismember(alldata_table_house_hem.roi,roi_str{roi_i}));
        alldata_table_house_hem_roi = alldata_table_house_hem(row_num,:);
        
        tableConditions.roi_test{roi_i} = unique(alldata_table_house_hem_roi.roi);
        tableConditions.condition_val = unique(alldata_table_house_hem_roi.condition);
        tableConditions.entryRoom_val = unique(alldata_table_house_hem_roi.entryRoom);
        tableConditions.room_cond_val = unique(alldata_table_house_hem_roi.room_cond);
        
        subj_list = unique(alldata_table_house_hem_roi.subj);
        
        all_z1_subj = []; all_z2_subj = [];
        all_z1_subj_v2 = []; all_z2_subj_v2 = [];
        for ii_s = 1:length(subj_list)
            row_num = find(ismember(alldata_table_house_hem_roi.subj,subj_list{ii_s}) &...
                ismember(alldata_table_house_hem_roi.room_cond,tableConditions.room_cond_val{2}) & ...
                ismember(alldata_table_house_hem_roi.condition,tableConditions.condition_val{2}));
            
            Z{1} = alldata_table_house_hem_roi.z(row_num);
            row_num = find(ismember(alldata_table_house_hem_roi.subj,subj_list{ii_s}) &...
                ismember(alldata_table_house_hem_roi.room_cond,tableConditions.room_cond_val{1}) & ...
                ismember(alldata_table_house_hem_roi.condition,tableConditions.condition_val{1}));
            
            all_z1_subj = [all_z1_subj Z{1}'];
            Z{2} = alldata_table_house_hem_roi.z(row_num);
            all_z2_subj = [all_z2_subj Z{2}'];
            pval_subj(ii_s,hem_i,roi_i) = ranksum(Z{1},Z{2});
            RSA_score_SUBJ_ROI(ii_s,hem_i,roi_i) = nanmedian(Z{1})-nanmedian(Z{2});
            
            
            % Second, less stringent RSA score - differences between houses
            % are sufficient?
            % SH/SR
            clear Z
            row_num = find(ismember(alldata_table_house_hem_roi.subj,subj_list{ii_s}) &...
                ismember(alldata_table_house_hem_roi.room_cond,tableConditions.room_cond_val{2}) & ...
                ismember(alldata_table_house_hem_roi.condition,tableConditions.condition_val{2}));
            % DH/SR
            Z{1} = alldata_table_house_hem_roi.z(row_num);
            row_num = find(ismember(alldata_table_house_hem_roi.subj,subj_list{ii_s}) &...
                ismember(alldata_table_house_hem_roi.room_cond,tableConditions.room_cond_val{2}) & ...
                ismember(alldata_table_house_hem_roi.condition,tableConditions.condition_val{1}));
            
            all_z1_subj_v2 = [all_z1_subj_v2 Z{1}'];
            Z{2} = alldata_table_house_hem_roi.z(row_num);
            all_z2_subj_v2 = [all_z2_subj_v2 Z{2}'];
            pval_subj_v2(ii_s,hem_i,roi_i) = ranksum(Z{1},Z{2});
            RSA_score_SR_SUBJ_ROI(ii_s,hem_i,roi_i) = nanmedian(Z{1})-nanmedian(Z{2});
            
        end
        pval_all_subj(hem_i,roi_i) = ranksum(all_z1_subj, all_z2_subj);
        pval_all_subj_v2(hem_i,roi_i) = ranksum(all_z1_subj_v2, all_z2_subj_v2);
        
        [R_space(hem_i,roi_i), corr_p_space(hem_i,roi_i)] = corr(RSA_score_SUBJ_ROI(:,hem_i,roi_i),HC_RC_perc_all');
        
        [R_rhit(hem_i,roi_i), corr_p_rhit(hem_i,roi_i)] = corr(RSA_score_SUBJ_ROI(:,hem_i,roi_i),rhit_perc_all');
        
        % Adding another less extreme condition SH/DH but SR
        [R_rhit_SR(hem_i,roi_i), corr_p_rhit_SR(hem_i,roi_i)] = corr(RSA_score_SR_SUBJ_ROI(:,hem_i,roi_i),rhit_perc_all');
        
        
    end
end

