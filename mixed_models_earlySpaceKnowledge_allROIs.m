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

row_num1 = find(ismember(alldata_table_house.entryLoc,'nearFrontDoor'));
alldata_table_house.border(row_num1) = {'border'};
row_num2 = find(ismember(alldata_table_house.passageLoc,'nearPassageDoor'));
alldata_table_house.border(row_num2) = {'border'};
rest_rows = ~ismember(1:length(alldata_table_house.passageLoc),[row_num1',row_num2']);
alldata_table_house.border(rest_rows) = {'other'};
unique(alldata_table_house.border)

row_num1 = find(ismember(alldata_table_house.entryRoom,'entryRoom'));
alldata_table_house.room_cond(row_num1) = {'sameRoom'};
row_num2 = find(ismember(alldata_table_house.entryRoom,'innerRoom'));
alldata_table_house.room_cond(row_num2) = {'sameRoom'};
rest_rows = ~ismember(1:length(alldata_table_house.room_cond),[row_num1',row_num2']);
alldata_table_house.room_cond(rest_rows) = {'diffRoom'};
unique(alldata_table_house.room_cond)

row_num1 = find(ismember(alldata_table_house.roi,'CA1_body'));
row_num2 = find(ismember(alldata_table_house.roi,'CA2_3_DG_body'));
alldata_table_house_ca1_ca23dg = alldata_table_house([row_num1(:)', row_num2(:)'],:) ;

formula = 'z~room_cond*roi + room_cond*hemi + roi*hemi + (1|subj)';
lme1 = fitlme(alldata_table_house_ca1_ca23dg,formula);
formula = 'z~room_cond*roi*hemi + (1|subj)';
lme2 = fitlme(alldata_table_house_ca1_ca23dg,formula);
results = compare(lme1,lme2);
disp('3-way interaction - entryLoc, roi, hem')
disp(results.pValue)

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

formula = 'z~entryLoc*roi + entryLoc*condition + condition*roi + (1|subj)';
lme1 = fitlme(alldata_table_house_left,formula);
formula = 'z~condition*roi*border + (1|subj)';
lme2 = fitlme(alldata_table_house_left,formula);
results = compare(lme1,lme2);
disp('3-way interaction - border, condition, roi')
disp(results.pValue)

formula = 'z~entryLoc*roi + entryLoc*condition + condition*roi + (1|subj)';
lme1 = fitlme(alldata_table_house_right,formula);
formula = 'z~condition*roi*entryLoc + (1|subj)';
lme2 = fitlme(alldata_table_house_right,formula);
results = compare(lme1,lme2);
disp('3-way interaction - main-door, condition, roi')
disp(results.pValue)

formula = 'z~entryLoc + roi + (1|subj)';
lme1 = fitlme(alldata_table_house_left,formula);
formula = 'z~roi*entryLoc + (1|subj)';
lme2 = fitlme(alldata_table_house_left,formula);
results = compare(lme1,lme2);
disp('LH - 2-way interaction - main-door, roi')
disp(results.pValue)

formula = 'z~entryLoc + roi + (1|subj)';
lme1 = fitlme(alldata_table_house_right,formula);
formula = 'z~roi*entryLoc + (1|subj)';
lme2 = fitlme(alldata_table_house_right,formula);
results = compare(lme1,lme2);
disp('RH - 2-way interaction - main-door, roi')
disp(results.pValue)

unique(alldata_table_house_left.hemi)
unique(alldata_table_house_left.roi)
unique(alldata_table_house_left.border)
unique(alldata_table_house_left.condition)
unique(alldata_table_house_left.room_cond)

%%
% All whithin the same house 
% Check various space-conditions 
pvalue_vec = [];  test_str = [];
for hem_i = 1:2
    for roi_i = 1:7
        if hem_i == 1
            alldata_table_house_hem = alldata_table_house_left;
        else
            alldata_table_house_hem = alldata_table_house_right;
        end
        
        row_num = find(ismember(alldata_table_house_hem.roi,roi_str{roi_i}));
        alldata_table_house_hem_roi = alldata_table_house_hem(row_num,:);
        
        row_num = find(ismember(alldata_table_house_hem_roi.condition,'differentVideo_sameHouse'));
        alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
        
        roi_test{roi_i} = unique(alldata_table_sameHouse_hem_roi.roi);
        condition_val = unique(alldata_table_sameHouse_hem_roi.condition);
        border_val = unique(alldata_table_sameHouse_hem_roi.border);
        entryRoom_val = unique(alldata_table_sameHouse_hem_roi.entryRoom);
        room_cond_val = unique(alldata_table_sameHouse_hem_roi.room_cond);
        
        test_i = 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_sameHouse_hem_roi,formula);
        formula = 'z~ entryRoom + (1|subj)';
        lme2 = fitlme(alldata_table_sameHouse_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'DV/SH R1/R2/DR';
        pvalue_vec(hem_i, roi_i,test_i) = results.pValue(2);
    
        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_sameHouse_hem_roi,formula);
        formula = 'z~ room_cond + (1|subj)';
        lme2 = fitlme(alldata_table_sameHouse_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'DV/SH SR/DR';
        pvalue_vec(hem_i, roi_i,test_i) = results.pValue(2);
        
        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_sameHouse_hem_roi,formula);
        formula = 'z~ entryLoc + (1|subj)';
        lme2 = fitlme(alldata_table_sameHouse_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'DV/SH  frontDoor';
        pvalue_vec(hem_i, roi_i,test_i) = results.pValue(2);
        
    end
end
test_str_DV_SH = test_str;
pvalue_vec_DV_SH = pvalue_vec;

%% Different videos - same/diff houses
% Check various space-conditions
pvalue_vec = [];  test_str = []; pvalue_ranksum_vec = [];
for hem_i = 1:3
    for roi_i = 1:7
        if hem_i == 1
            alldata_table_house_hem = alldata_table_house_left;
        elseif hem_i == 2
            alldata_table_house_hem = alldata_table_house_right;
        else
             alldata_table_house_hem = alldata_table_house;
        end
        row_num = find(ismember(alldata_table_house_hem.roi,roi_str{roi_i}));
        alldata_table_house_hem_roi = alldata_table_house_hem(row_num,:);
        
        tableConditions.roi_test{roi_i} = unique(alldata_table_house_hem_roi.roi);
        tableConditions.condition_val = unique(alldata_table_house_hem_roi.condition);
        tableConditions.border_val = unique(alldata_table_house_hem_roi.border);
        tableConditions.entryRoom_val = unique(alldata_table_house_hem_roi.entryRoom);
        tableConditions.room_cond_val = unique(alldata_table_house_hem_roi.room_cond);
        
        test_i = 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_hem_roi,formula);
        formula = 'z~ entryRoom + (1|subj)';
        lme2 = fitlme(alldata_table_house_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'DV/SH R1/R2/DR';
        pvalue_vec(hem_i,roi_i,test_i) = results.pValue(2);

                
        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_hem_roi,formula);
        formula = 'z~ condition + (1|subj)';
        lme2 = fitlme(alldata_table_house_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'pure space 1 - same vs diff house';
        pvalue_vec(hem_i, roi_i,test_i) = results.pValue(2);
        
        test_i = test_i + 1;
        formula = 'z~ room_cond + condition + (1|subj)';
        lme1 = fitlme(alldata_table_house_hem_roi,formula);
        formula = 'z~ condition*room_cond + (1|subj)';
        lme2 = fitlme(alldata_table_house_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'pure space 2 - interaction - house and room';
        pvalue_vec(hem_i, roi_i,test_i) = results.pValue(2);
        clear AAA
        for ii_r = 1:2
            rows = find(ismember(alldata_table_house_hem_roi.condition,tableConditions.condition_val{ii_r}) &...
                        ismember(alldata_table_house_hem_roi.room_cond,tableConditions.room_cond_val{ii_r}) );
            AAA{ii_r} = alldata_table_house_hem_roi.z(rows);
        end
        pvalue_ranksum_vec(hem_i,roi_i,test_i) = ranksum(AAA{1},AAA{2});

        
        row_num = find(ismember(alldata_table_house_hem_roi.entryLoc,'nearFrontDoor'));
        alldata_table_house_frontDoor_hem_roi = alldata_table_house_hem_roi(row_num,:);
        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_frontDoor_hem_roi,formula);
        formula = 'z~ condition + (1|subj)';
        lme2 = fitlme(alldata_table_house_frontDoor_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'doorway SH/DH';
        pvalue_vec(hem_i,roi_i,test_i) = results.pValue(2);
        clear AAA
        for ii_r = 1:2
            rows = find(ismember(alldata_table_house_frontDoor_hem_roi.condition,tableConditions.condition_val{ii_r}));
            AAA{ii_r} = alldata_table_house_frontDoor_hem_roi.z(rows);
        end
        pvalue_ranksum_vec(hem_i,roi_i,test_i) = ranksum(AAA{1},AAA{2});

        row_num = find(ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
        alldata_table_house_frontDoor_hem_roi = alldata_table_house_hem_roi(row_num,:);
        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_frontDoor_hem_roi,formula);
        formula = 'z~ condition + (1|subj)';
        lme2 = fitlme(alldata_table_house_frontDoor_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'entryRoom SH/DH';
        pvalue_vec(hem_i,roi_i,test_i) = results.pValue(2);
        clear AAA
        for ii_r = 1:2
            rows = find(ismember(alldata_table_house_frontDoor_hem_roi.condition,tableConditions.condition_val{ii_r}));
            AAA{ii_r} = alldata_table_house_frontDoor_hem_roi.z(rows);
        end
        pvalue_ranksum_vec(hem_i,roi_i,test_i) = ranksum(AAA{1},AAA{2});

        % Test items in entry room vs all others
        row_num = find(ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
        alldata_table_house_hem_roi.entryRoom_vs_all(row_num) = {'entryRoom'};
        row_num = find(~ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
        alldata_table_house_hem_roi.entryRoom_vs_all(row_num) = {'other'};
        entryRoom_vs_all_cond = unique(alldata_table_house_hem_roi.entryRoom_vs_all);
        
        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_hem_roi,formula);
        formula = 'z~ entryRoom_vs_all + (1|subj)';
        lme2 = fitlme(alldata_table_house_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'ALL entryRoom vs other';
        pvalue_vec(hem_i,roi_i,test_i) = results.pValue(2);
        clear AAA
        for ii_r = 1:2
            rows = find(ismember(alldata_table_house_hem_roi.entryRoom_vs_all,entryRoom_vs_all_cond{ii_r}));
            AAA{ii_r} = alldata_table_house_hem_roi.z(rows);
        end
        pvalue_ranksum_vec(hem_i,roi_i,test_i) = ranksum(AAA{1},AAA{2});
        

        row_num = find(ismember(alldata_table_house_hem_roi.border,'border'));
        alldata_table_house_frontDoor_hem_roi = alldata_table_house_hem_roi(row_num,:);
        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_frontDoor_hem_roi,formula);
        formula = 'z~ condition + (1|subj)';
        lme2 = fitlme(alldata_table_house_frontDoor_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'border SH/DH';
        pvalue_vec(hem_i,roi_i,test_i) = results.pValue(2);
        clear AAA
        for ii_r = 1:2
            rows = find(ismember(alldata_table_house_frontDoor_hem_roi.condition,tableConditions.condition_val{ii_r}));
            AAA{ii_r} = alldata_table_house_frontDoor_hem_roi.z(rows);
        end
        pvalue_ranksum_vec(hem_i,roi_i,test_i) = ranksum(AAA{1},AAA{2});

        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_hem_roi,formula);
        formula = 'z~ room_cond + (1|subj)';
        lme2 = fitlme(alldata_table_house_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'ALL SR/DR';
        pvalue_vec(hem_i,roi_i,test_i) = results.pValue(2);
        clear AAA
        for ii_r = 1:2
            rows = find(ismember(alldata_table_house_frontDoor_hem_roi.condition,tableConditions.condition_val{ii_r}));
            AAA{ii_r} = alldata_table_house_frontDoor_hem_roi.z(rows);
        end
        pvalue_ranksum_vec(hem_i,roi_i,test_i) = ranksum(AAA{1},AAA{2});
        
        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_hem_roi,formula);
        formula = 'z~ entryLoc + (1|subj)';
        lme2 = fitlme(alldata_table_house_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'ALL frontdoor/other';
        pvalue_vec(hem_i,roi_i,test_i) = results.pValue(2);

        test_i = test_i + 1;
        formula = 'z~ 1  + (1|subj)';
        lme1 = fitlme(alldata_table_house_hem_roi,formula);
        formula = 'z~ border + (1|subj)';
        lme2 = fitlme(alldata_table_house_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'ALL any-doorway/other';
        pvalue_vec(hem_i,roi_i,test_i) = results.pValue(2);
        
    end
end
test_str_DV = test_str;
pvalue_vec_DV = pvalue_vec;
stats_dir = 'E:\Dropbox\RanganathLab\spatcon_2\stats\';
save(fullfile(stats_dir,'mixed_effects_models'),'test_str_DV','pvalue_vec_DV','pvalue_ranksum_vec')

%%
if PLOT_fig

    for hem_i = 1:2
        
        if hem_i == 1
            title_str1 = sprintf('LH - entryRoom_roi');
            title_str2 = sprintf('LH - pureSpace_roi');
        else
            title_str1 = sprintf('RH - entryRoom_roi');
            title_str2 = sprintf('RH - pureSpace_roi');
            
        end
        f0 = newA4figure(title_str1);
        f1 = newA4figure(title_str2);
        
        roi_list = [1,2,4,5:7];
        for roi_ind = 1:length(roi_list);
            
            roi_i = roi_list(roi_ind);
            
            if hem_i == 1
                alldata_table_house_hem = alldata_table_house_left;
            else
                alldata_table_house_hem = alldata_table_house_right;
            end
            row_num = find(ismember(alldata_table_house_hem.roi,roi_str{roi_i}));
            alldata_table_house_hem_roi = alldata_table_house_hem(row_num,:);
            
            % entryRoom
            rows = find( ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
            Z{1,1} =  alldata_table_house_hem_roi.z(rows);
            rows = find(ismember(alldata_table_house_hem_roi.entryRoom,'innerRoom'));
            Z{1,2}=  alldata_table_house_hem_roi.z(rows);
            rows = find(ismember(alldata_table_house_hem_roi.entryRoom,'otherLoc') );
            Z{2,1} =  alldata_table_house_hem_roi.z(rows);
            rows = find(ismember(alldata_table_house_hem_roi.entryRoom,'otherLoc') );
            Z{2,2} =  alldata_table_house_hem_roi.z(rows);
            
            figure(f0); subplot(length(roi_list),2,(roi_ind-1)*2+1)
            xlabel_str = {'SR - 1','DR','SR - 2','DR'};
            plot_4by4_overlap_distributions(Z,xlabel_str)
            if hem_i == 2
                title(sprintf('RH - %s',regexprep(roi_str{roi_i},'_',' ')))
            else
                title(sprintf('LH - %s',regexprep(roi_str{roi_i},'_',' ')))
            end
            
            clear Z
            
            clear Z
            rows = find( ismember(alldata_table_house_hem_roi.entryLoc,'nearFrontDoor'));
            Z{1} =  alldata_table_house_hem_roi.z(rows);
            rows = find(ismember(alldata_table_house_hem_roi.entryLoc,'otherLoc'));
            Z{2}=  alldata_table_house_hem_roi.z(rows);
            
            figure(f0); subplot(length(roi_list),2,(roi_ind-1)*2+2)
            xlabel_str = {'EntryDoor','other'};
            plot_2by2_distributions(Z,xlabel_str)
            if hem_i == 2
                title(sprintf('RH - %s',regexprep(roi_str{roi_i},'_',' ')))
            else
                title(sprintf('LH - %s',regexprep(roi_str{roi_i},'_',' ')))
            end
            axis([0,3,0,inf])
            
            
            %% Figure 2
            % same house same room
            clear Z
            row_num = find(ismember(alldata_table_house_hem_roi.condition,'differentVideo_sameHouse'));
            alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            rows = find( ismember(alldata_table_sameHouse_hem_roi.room_cond,'sameRoom'));
            Z{1} =  alldata_table_sameHouse_hem_roi.z(rows);
            
            row_num = find(ismember(alldata_table_house_hem_roi.condition,'differentHouse_eitherRoom'));
            alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            rows = find( ismember(alldata_table_sameHouse_hem_roi.room_cond,'diffRoom'));
            Z{2} =  alldata_table_sameHouse_hem_roi.z(rows);
            
            figure(f1); subplot(length(roi_list),2,(roi_ind-1)*2+1)
            xlabel_str = {'SH/SR','DH/DR'};
            plot_2by2_distributions(Z,xlabel_str)
            if hem_i == 2
                title(sprintf('RH - %s',regexprep(roi_str{roi_i},'_',' ')))
            else
                title(sprintf('LH - %s',regexprep(roi_str{roi_i},'_',' ')))
            end
            [pval_space(hem_i, roi_i),h] = ranksum(Z{1},Z{2});
            axis([0,3,0,inf])
            
            % focus on doorway area
            clear Z
            row_num = find(ismember(alldata_table_house_hem_roi.condition,'differentVideo_sameHouse'));
            alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            rows = find( ismember(alldata_table_sameHouse_hem_roi.entryRoom,'entryRoom'));
            Z{1} =  alldata_table_sameHouse_hem_roi.z(rows);
            
            row_num = find(ismember(alldata_table_house_hem_roi.condition,'differentHouse_eitherRoom'));
            alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            rows = find( ismember(alldata_table_sameHouse_hem_roi.entryRoom,'entryRoom'));
            Z{2} =  alldata_table_sameHouse_hem_roi.z(rows);
            
            figure(f1); subplot(length(roi_list),2,(roi_ind-1)*2+2)
            xlabel_str = {'SH/R1','DH/R1'};
            plot_2by2_distributions(Z,xlabel_str)
            if hem_i == 2
                title(sprintf('RH - %s',regexprep(roi_str{roi_i},'_',' ')))
            else
                title(sprintf('LH - %s',regexprep(roi_str{roi_i},'_',' ')))
            end
            [pval_doorway(hem_i, roi_i),h] = ranksum(Z{1},Z{2});
            axis([0,3,0,inf])
            
        end
        
    end
end

%%
if PLOT_fig2

    for hem_i = 1:2
        
        if hem_i == 1
            title_str1 = sprintf('LH - entryRoom_roi');
            title_str2 = sprintf('LH - pureSpace_roi');
        else
            title_str1 = sprintf('RH - entryRoom_roi');
            title_str2 = sprintf('RH - pureSpace_roi');
            
        end
        f0 = newA4figure(title_str1);
        f1 = newA4figure(title_str2);
        
        roi_list = [1,2,4,5:7];
        for roi_ind = 1:length(roi_list);
            
            roi_i = roi_list(roi_ind);
            
            if hem_i == 1
                alldata_table_house_hem = alldata_table_house_left;
            else
                alldata_table_house_hem = alldata_table_house_right;
            end
            row_num = find(ismember(alldata_table_house_hem.roi,roi_str{roi_i}));
            alldata_table_house_hem_roi = alldata_table_house_hem(row_num,:);
            
            % entryRoom
            clear Z
            rows = find( ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
            Z{1,1} =  alldata_table_house_hem_roi.z(rows);
            rows = find( ~ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
            Z{1,2}=  alldata_table_house_hem_roi.z(rows);
            
            figure(f0); subplot(length(roi_list),2,(roi_ind-1)*2+1)
            xlabel_str = {'entryRoom','other'};
            plot_2by2_distributions(Z,xlabel_str)
          
            if hem_i == 2
                title(sprintf('RH - %s',regexprep(roi_str{roi_i},'_',' ')))
            else
                title(sprintf('LH - %s',regexprep(roi_str{roi_i},'_',' ')))
            end
            axis([0,3,0,inf])
            
            clear Z
            rows = find( ismember(alldata_table_house_hem_roi.entryLoc,'nearFrontDoor'));
            Z{1} =  alldata_table_house_hem_roi.z(rows);
            rows = find(ismember(alldata_table_house_hem_roi.entryLoc,'otherLoc'));
            Z{2}=  alldata_table_house_hem_roi.z(rows);
            
            figure(f0); subplot(length(roi_list),2,(roi_ind-1)*2+2)
            xlabel_str = {'EntryDoor','other'};
            plot_2by2_distributions(Z,xlabel_str)
            if hem_i == 2
                title(sprintf('RH - %s',regexprep(roi_str{roi_i},'_',' ')))
            else
                title(sprintf('LH - %s',regexprep(roi_str{roi_i},'_',' ')))
            end
            axis([0,3,0,inf])
            
            
            % Figure 2
            % same house same room
            clear Z
            row_num = find(ismember(alldata_table_house_hem_roi.condition,'differentVideo_sameHouse') &...
                            ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
            alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            Z{1} =  alldata_table_sameHouse_hem_roi.z;

            row_num = find(~ismember(alldata_table_house_hem_roi.condition,'differentVideo_sameHouse') |...
                           ismember(alldata_table_house_hem_roi.entryRoom,'otherLoc'));
            alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            Z{2} =  alldata_table_sameHouse_hem_roi.z;
            
            figure(f1); subplot(length(roi_list),2,(roi_ind-1)*2+1)
            xlabel_str = {'SH/entryRoom','DH/DR'};
            plot_2by2_distributions(Z,xlabel_str)
            if hem_i == 2
                title(sprintf('RH - %s',regexprep(roi_str{roi_i},'_',' ')))
            else
                title(sprintf('LH - %s',regexprep(roi_str{roi_i},'_',' ')))
            end
            [pval_space(hem_i, roi_i),h] = ranksum(Z{1},Z{2});
            axis([0,3,0,inf])
            
            % focus on doorway area
            clear Z
            row_num = find(ismember(alldata_table_house_hem_roi.condition,'differentVideo_sameHouse'));
            alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            rows = find( ismember(alldata_table_sameHouse_hem_roi.entryLoc,'nearFrontDoor'));
            Z{1} =  alldata_table_sameHouse_hem_roi.z(rows);
            
            row_num = find(ismember(alldata_table_house_hem_roi.condition,'differentHouse_eitherRoom'));
            alldata_table_sameHouse_hem_roi = alldata_table_house_hem_roi(row_num,:);
            rows = find( ismember(alldata_table_sameHouse_hem_roi.entryLoc,'nearFrontDoor'));
            Z{2} =  alldata_table_sameHouse_hem_roi.z(rows);
            
            figure(f1); subplot(length(roi_list),2,(roi_ind-1)*2+2)
            xlabel_str = {'SH/front-door','DH/front-door'};
            plot_2by2_distributions(Z,xlabel_str)
            if hem_i == 2
                title(sprintf('RH - %s',regexprep(roi_str{roi_i},'_',' ')))
            else
                title(sprintf('LH - %s',regexprep(roi_str{roi_i},'_',' ')))
            end
            [pval_doorway(hem_i, roi_i),h] = ranksum(Z{1},Z{2});
            axis([0,3,0,inf])
            
        end
        
    end
end


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

% Hypothesis - higher sensitivity for space will be correlated with better
% judgement, 'memorability'.

% define an RSA grade per subject for the difference between sameH/sameR
% and diff house items.

pvalue_vec = [];  test_str = []; pvalue_ranksum_vec = [];
clear R_space R_rhit corr_p_space corr_p_rhit
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
        tableConditions.border_val = unique(alldata_table_house_hem_roi.border);
        tableConditions.entryRoom_val = unique(alldata_table_house_hem_roi.entryRoom);
        tableConditions.room_cond_val = unique(alldata_table_house_hem_roi.room_cond);
        
        subj_list = unique(alldata_table_house_hem_roi.subj);
        
        for ii_s = 1:length(subj_list)
            row_num = find(ismember(alldata_table_house_hem_roi.subj,subj_list{ii_s}) &...
                ismember(alldata_table_house_hem_roi.room_cond,tableConditions.room_cond_val{2}) & ...
                ismember(alldata_table_house_hem_roi.condition,tableConditions.condition_val{2}));
            
            Z{1} = alldata_table_house_hem_roi.z(row_num);
            row_num = find(ismember(alldata_table_house_hem_roi.subj,subj_list{ii_s}) &...
                ismember(alldata_table_house_hem_roi.room_cond,tableConditions.room_cond_val{1}) & ...
                ismember(alldata_table_house_hem_roi.condition,tableConditions.condition_val{1}));
            Z{2} = alldata_table_house_hem_roi.z(row_num);
            pval_subj(ii_s,hem_i,roi_i) = ranksum(Z{1},Z{2});
            RSA_score_SUBJ_ROI(ii_s,hem_i,roi_i) = nanmedian(Z{1})-nanmedian(Z{2});
        end
        [R_space(hem_i,roi_i), corr_p_space(hem_i,roi_i)] = corr(RSA_score_SUBJ_ROI(:,hem_i,roi_i),HC_RC_perc_all');

        [R_rhit(hem_i,roi_i), corr_p_rhit(hem_i,roi_i)] = corr(RSA_score_SUBJ_ROI(:,hem_i,roi_i),rhit_perc_all');
    end
end


dataFolder = 'E:\Dropbox\RanganathLab\spatcon_2\figures\dataForFigure';
save(fullfile(dataFolder, 'RSA_score_SUBJ_ROI'), 'RSA_score_SUBJ_ROI','HC_RC_perc_all','roi_str');

%% Relationship to behav - 
newA4figure('RSA_diff_corr_behav_space')
ccmap = brewermap(length(HC_RC_perc_all)+5,'OrRd');

for ii_plot = 1:3
ax1 = subplot(1,3,ii_plot);

if ii_plot == 1
    hem_i = 2;
    roi_i = 4;
elseif ii_plot == 2
    hem_i = 2;
    roi_i = 7;
else
    hem_i = 2;
    roi_i = 6;
end

[B,ind] = sort(HC_RC_perc_all);
C = RSA_score_SUBJ_ROI(ind,hem_i,roi_i);

hold all
for ii = 1:length(HC_RC_perc_all)
    plot(C(ii),100*B(ii),'o','color',ccmap(5+ii,:),'linewidth',2)
    if pval_subj(ii,hem_i,roi_i)<0.05
        plot(C(ii),100*B(ii),'k*','markersize',5)
    end
end
axis square
set(gca,'fontsize',12)
set(gca,'ylim',100*[0 1])
set(gca,'YTick',100*[0 0.5 1])
set(gca,'Xlim',[-20  5]*10-3)
set(gca,'Xlim',[-20  20]*10^-3)
set(gca,'XTick',[-20 0 20]*10^-3)
hold all
plot([0 0],get(gca,'ylim'),'--k')
XLIM = get(gca,'xlim');
YLIM = get(gca,'ylim');
text(0.5*diff(XLIM)+XLIM(1),YLIM(2)+0.03*diff(YLIM),...
   sprintf('%s, R = %2.2f, p = %2.2e',strrep(roi_str{roi_i},'_',' '), R_space(hem_i,roi_i),corr_p_space(hem_i,roi_i)),...
   'HorizontalAlignment','center');
xlabel('RSA spatial sensitivity score')
ylabel('Correct placement in post-scan test (%)')

coefs = polyfit(C, 100*B', 1);
% Plot line
ah = refline(coefs); 
ah.Color = 'k';
set(gca,'ylim',100*[0 1])

end


newA4figure('RSA_diff_corr_behav_rhit')
ccmap = brewermap(length(rhit_perc_all)+5,'OrRd');

for ii_plot = 1:3
ax1 = subplot(1,3,ii_plot);

if ii_plot == 1
    hem_i = 2;
    roi_i = 4;
elseif ii_plot == 2
    hem_i = 2;
    roi_i = 7;
else
    hem_i = 2;
    roi_i = 6;
end

[B,ind] = sort(rhit_perc_all);
C = RSA_score_SUBJ_ROI(ind,hem_i,roi_i);

hold all
for ii = 1:length(rhit_perc_all)
    plot(C(ii),100*B(ii),'o','color',ccmap(5+ii,:),'linewidth',2)
    if pval_subj(ii,hem_i,roi_i)<0.05
        plot(C(ii),100*B(ii),'k*','markersize',5)
    end
end
axis square
set(gca,'fontsize',12)
set(gca,'ylim',100*[0 1])
set(gca,'YTick',100*[0 0.5 1])
set(gca,'Xlim',[-20  5]*10-3)
set(gca,'Xlim',[-20  20]*10^-3)
set(gca,'XTick',[-20 0 20]*10^-3)
hold all
plot([0 0],get(gca,'ylim'),'--k')
XLIM = get(gca,'xlim');
YLIM = get(gca,'ylim');
text(0.5*diff(XLIM)+XLIM(1),YLIM(2)+0.03*diff(YLIM),...
   sprintf('%s, R = %2.2f, p = %2.2e',strrep(roi_str{roi_i},'_',' '), R_rhit(hem_i,roi_i),corr_p_rhit(hem_i,roi_i)),...
   'HorizontalAlignment','center');
xlabel('RSA spatial sensitivity score')
ylabel('Recollection (%)')

coefs = polyfit(C, 100*B', 1);
% Plot line
ah = refline(coefs); 
ah.Color = 'k';
set(gca,'ylim',100*[0 1])

end

PrintActiveFigs('E:\Dropbox\RanganathLab\spatcon_2\figures\LINKS');

