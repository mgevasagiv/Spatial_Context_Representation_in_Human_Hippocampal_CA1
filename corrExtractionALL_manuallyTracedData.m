% Generate correlation matrices that would enter the Matlab pipeline (based on original code that generated these matrices in python)
% This is similar to corrExtractionALL_publishedData.m, this code uploades additional MTL areas (outside
% the hippocampal subfields and EC)
clear all; 


%% Loading
% subject sepcific onset files are extracted using stimuliExtraction.m
onsetFiles_root = 'E:\spatcon_hdz\data_p\onsetFiles\';
analyzed_mri_dir = 'E:\spatcon_hdz\data_p\analyzed_mri_try1/'; % based on OSF downloads

% column numbers -
SUBJ = 1;
ROI = 2;
HEM = 3;
COND = 4;
ROW_NAME = 5;
COL_NAME = 6;
number_of_trials_th  = 20;

RIGHT = 1;
LEFT = 2;

SH_SR = 1;
DH_DR = 2;
SV_SR = 3;
DV_SR = 4;

roi_dir = 'manual';

subjects = 1:30;
exclude_subjects = [3:5,15:17];
subjects = subjects(~ismember(subjects,exclude_subjects));

%%
for subj_i = 1:length(subjects)
    
    subj = subjects(subj_i);
    disp(subj)
    
    for hem_i = 1:2
        
        if hem_i == 1
            hemi_label = 'left';
            rois_of_most_interest = { 'pmERC_L', 'alERC_L',...
                                        'PHC_Ant_L', 'PRC_L'};
        else
            hemi_label = 'right';
            rois_of_most_interest = {'pmERC_R','alERC_R', ...
                                         'PHC_Ant_R', 'PRC_R'};
        end
        
        clear pattern_corr_orig fisherz_pattern_corr pattern_corr_orig
        for roi_i = 1:length(rois_of_most_interest)
            
            roi_str = rois_of_most_interest{roi_i};
            
            % 
            fileName = fullfile(analyzed_mri_dir,sprintf('s%03d',subj),'ROIs',roi_dir,...
                sprintf('br%s_pattern_corr_no_outlier_trials_all_runs.mat',roi_str));
            fileNameIDs = fullfile(analyzed_mri_dir,sprintf('s%03d',subj),'ROIs',roi_dir,...
                sprintf('br%s_pattern_mtx_ids_no_outlier_trials_all_runs.mat',roi_str));
            
            if isempty(dir(fileName))
                disp([fileName,' is missing'])
                continue
            end
            
            
            mm = matfile(fileName);
            pattern_corr_orig = mm.pattern_corr;
            pattern_corr = pattern_corr_orig;
            pattern_corr(isnan(pattern_corr_orig)) = 0;
            % to reduce skewness(asymmetry derived from the definition of riin the interval [-1 1]), 
            % rvalues were transformed to zvalues via the Fisher procedure
            % (this is based on Cohen 1988, Statistical power analysis for the behavioral sciences)
            fisherz_pattern_corr = 0.5*log((1+pattern_corr)./(1-pattern_corr)); % fisher z-transform (to avoid high values)
            mm_id = matfile(fileNameIDs);
            pattern_ids_all_runs = mm_id.pattern_ids_all_runs;
                      
            folderName = fullfile(onsetFiles_root, sprintf('s%03d',subj));
            try
                stimuliFile = fullfile(folderName,sprintf('stimuliList_s%03d.mat',subj));
                mm = matfile(stimuliFile);
                objectID = mm.objectID;
                REF_ID = mm.REF_ID;
                
            catch
                warning('onset file not found')
            end
            
            
            if ((hem_i == 1) && (roi_i == 1)) || (subj == 22 & roi_i == 3)
                clear object_num mem1_str mem2_str 
                for ii_s = 1:length(pattern_ids_all_runs)
                    str = pattern_ids_all_runs{ii_s};
                    C = strsplit(str,'_');
                    if strcmp(C{5},'ExcludeTrial') || strcmp(C{5},'CR') || strcmp(C{5},'Miss') || strcmp(C{5},'FA')
                        object_num(ii_s) = NaN;
                        mem1_str{ii_s} = NaN;
                        mem2_str{ii_s} = NaN;
                        continue
                    end
                    run_i = str2num(C{2});
                    trial_i = str2num(C{4});
                    if strcmp(C{7},'Lure') || strcmp(C{8},'Lure')
                        house_i = 0;
                    elseif strcmp(C{5}(1),'B')
                        house_i = 1;
                    else
                        house_i = -1;
                    end
                    index = find(objectID(:,REF_ID.trialNum) == trial_i & objectID(:,REF_ID.house) == house_i & ...
                        objectID(:,REF_ID.BlockID) == run_i);
                    if length(index) > 1
                        error('not specific enough')
                    else
                        object_num(ii_s) = objectID(index,REF_ID.object_num);
                        mem1_str{ii_s} = C{6};
                        mem2_str{ii_s} = C{7};
                        
                    end
                end
                
                object_info_all_rois{subj}.object_num = object_num;
                object_info_all_rois{subj}.mem1_str = mem1_str;
                object_info_all_rois{subj}.mem2_str = mem2_str;
                object_info_all_rois{subj}.pattern_ids_all_runs = pattern_ids_all_runs;
            end
            
            pattern_corr_all_rois{subj,hem_i, roi_i} = pattern_corr_orig;
            fz_pattern_corr_all_rois{subj,hem_i, roi_i} = fisherz_pattern_corr;
            
        end
    end
end

save('E:\spatcon_hdz\data_p\corrAnalysis_publishedData\pattern_corr_manualTracing_zscored_all_subjects.mat','object_info_all_rois','pattern_corr_all_rois','rois_of_most_interest')

%%
% Plot raw data and zscored data per subject for all ROIs - note that the
% most prominent correlations are whithin run and these are removed (can
% they be used somehow?)

% Plot per subject
figureFolder1 = 'E:\spatcon_hdz\data_p\figures\correlationFigures\rawR';
figureFolder2 = 'E:\spatcon_hdz\data_p\figures\correlationFigures\fisherZ';
hem_i = 1;
rois_of_most_interest= { 'pmERC_L', 'alERC_L',...
    'PHC_Ant_L', 'PRC_L'};
     
for subj_i = 1:length(object_info_all_rois)
    if isempty(object_info_all_rois{subj_i})
        continue
    end
    
    disp(subj_i)

    title_f = sprintf('s%03d_rawR_leftH_all_manualROIs',subj_i);
    f0 = newA4figure(title_f);
    for roi_i = 1:length(rois_of_most_interest)
        
           title_s = sprintf('%s',rois_of_most_interest{roi_i});
            
           subplot(2,3,roi_i)
           MAT = pattern_corr_all_rois{subj_i,hem_i, roi_i};
           imagesc(MAT)
           title(['roi: ',strrep(title_s,'_',' ')])
           colorbar
    end
    saveas(f0,fullfile(figureFolder1,[title_f,'.jpg']))
    close(f0)
    
    
    title_f = sprintf('s%03d_fisherZ_R_leftH_all_manualROIs',subj_i);
    f0 = newA4figure(title_f);
    for roi_i = 1:length(rois_of_most_interest)
        
           title_s = sprintf('%s',rois_of_most_interest{roi_i});
            
           subplot(2,3,roi_i)
           MAT = fz_pattern_corr_all_rois{subj_i,hem_i, roi_i};
           imagesc(MAT)
           title(['roi: ',strrep(title_s,'_',' ')])
           colorbar
    end
    saveas(f0,fullfile(figureFolder2,[title_f,'.jpg']))
    close(f0)
    
end

