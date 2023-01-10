function gen_spatial_RSA_fig()

% prep data
mixed_models_add_behav_for_fig()

condition_val = unique(alldata_table_house.condition);
disp(condition_val{2})

% STATS 1 for MS - 
% restrict for trials from different videos exploring the same house -
row_num = find(ismember(alldata_table_house.condition,condition_val{2}));
alldata_table_SameHouse = alldata_table_house(row_num,:);
    
formula = 'z~ roi+entryRoom +(1|subj)';
lme1 = fitlme(alldata_table_SameHouse,formula);
formula = 'z~entryRoom*roi + (1|subj)';
lme2 = fitlme(alldata_table_SameHouse,formula);
results = compare(lme1,lme2);
disp('Same house (diff video) - 2-way interaction (both hem)- roi, room') 
disp(results.pValue) % P = 2.628e-05

formula = 'z~ entryRoom + hemi+(1|subj)';
lme1 = fitlme(alldata_table_house,formula);
formula = 'z~entryRoom*hemi + (1|subj)';
lme2 = fitlme(alldata_table_house,formula);
results = compare(lme1,lme2);
disp('2-way interaction - hemi,room') % P = 0.64353
disp(results.pValue)

% formula = 'z~ entryRoom + hemi+(1|subj)';
% lme1 = fitlme(alldata_table_house_ca1_ca23dg,formula);
% formula = 'z~entryRoom*hemi + (1|subj)';
% lme2 = fitlme(alldata_table_house_ca1_ca23dg,formula);
% results = compare(lme1,lme2);
% disp('2-way interaction - hemi,room') % P = 0.4
% disp(results.pValue)

formula = 'z~ entryRoom*hemi + roi*entryRoom + hemi*roi+(1|subj)';
lme1 = fitlme(alldata_table_house,formula);
formula = 'z~entryRoom*roi*hemi + (1|subj)';
lme2 = fitlme(alldata_table_house,formula);
results = compare(lme1,lme2);
disp('3-way interaction - hemi, roi, room') % P = 0.0025075
disp(results.pValue)


formula = 'z~ entryRoom*hemi + roi*entryRoom + hemi*roi+(1|subj)';
lme1 = fitlme(alldata_table_house_ca1_ca23dg,formula);
formula = 'z~entryRoom*roi*hemi + (1|subj)';
lme2 = fitlme(alldata_table_house_ca1_ca23dg,formula);
results = compare(lme1,lme2);
disp('3-way interaction - hemi, roi, room') % P = 0.43
disp(results.pValue)

formula = 'z~ roi+entryRoom +(1|subj)';
lme1 = fitlme(alldata_table_house_ca1_ca23dg,formula);
formula = 'z~entryRoom*roi + (1|subj)';
lme2 = fitlme(alldata_table_house_ca1_ca23dg,formula);
results = compare(lme1,lme2);
disp('2-way interaction (both hem)- roi, room') 
disp(results.pValue) % P = 6.9882e-06

formula = 'z~ roi+entryRoom +(1|subj)';
lme1 = fitlme(alldata_table_house_left,formula);
formula = 'z~entryRoom*roi + (1|subj)';
lme2 = fitlme(alldata_table_house_left,formula);
results = compare(lme1,lme2);
disp('2-way interaction (left hem)- roi, room') % 0.00063446
disp(results.pValue)

formula = 'z~ roi+entryRoom +(1|subj)';
lme1 = fitlme(alldata_table_house_right,formula);
formula = 'z~entryRoom*roi + (1|subj)';
lme2 = fitlme(alldata_table_house_right,formula);
results = compare(lme1,lme2);
disp('2-way interaction (right hem)- roi, room') % 8.5017e-06
disp(results.pValue)

% Just space
formula = 'z~ 1 +(1|subj)';
lme1 = fitlme(alldata_table_house_ca1_ca23dg,formula);
formula = 'z~ spaceCond_joint + (1|subj)';
lme2 = fitlme(alldata_table_house_ca1_ca23dg,formula);
results = compare(lme1,lme2);
disp('just space')
disp(results.pValue) % n.s.



%%
% Nice looking figures for Manuscript figures
x_w = 0.08;
pos(1,:) = [0.1, 0.1, x_w, 0.15];
pos(2,:) = [0.27, 0.1, x_w, 0.15];
pos(3,:) = [0.44, 0.1, x_w, 0.15];
pos(4,:) = [0.7, 0.1, x_w, 0.15];

cmap = brewermap(5,'set1');
barcolor(1,:) = cmap(5,:);
barcolor(2,:) = cmap(3,:);
barcolor(3,:) = [0.6 0.6 0.6];


title_str1 = 'fig2b';
f0 = newA4figure(title_str1);
pval_room = [];
roi_list = [1,2];
xlabel_str = {'entryRoom','innerRoom','cross room'};


for roi_ind = 1:length(roi_list)
    
    roi_i = roi_list(roi_ind);
    
    % restrict to PS of analyzed ROI - 
    row_num = find(ismember(alldata_table_house.roi,roi_str{roi_i}));
    alldata_table_house_roi = alldata_table_house(row_num,:);

    % restrict for trials from different videos exploring the same house - 
    row_num = find(ismember(alldata_table_house_roi.condition,condition_val{2}));
    alldata_table_SameHouse_roi = alldata_table_house_roi(row_num,:);
    
    if roi_ind == 1
        YTICKS = [0 4]*10^-3;
    elseif roi_ind == 2
        YTICKS = [0 8]*10^-3;
    elseif roi_ind == 3
        YTICKS = [0 ,3.5]*10^-2;
    elseif roi_ind == 4
        YTICKS = [0 ,10]*10^-3;
    end
    
     
    % main panel - same house entryRoom
    clear Z
    rows = find( ismember(alldata_table_SameHouse_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_SameHouse_roi.z(rows);
    rows = find( ismember(alldata_table_SameHouse_roi.entryRoom,'innerRoom'));
    Z{1,2}=  alldata_table_SameHouse_roi.z(rows);
    rows = find( ismember(alldata_table_SameHouse_roi.entryRoom,'otherLoc'));
    Z{1,3}=  alldata_table_SameHouse_roi.z(rows);
  
    y = [Z{1,1}',Z{1,2}',Z{1,3}'];
    g = [1*ones(1,length(Z{1,1})),2*ones(1,length(Z{1,2})),3*ones(1,length(Z{1,3}))];
    [p_anovan(roi_i),tbl,stats] = anovan(y,{g});
    figure; multcompare(stats)
    
    formula = 'z~ 1 +(1|subj)';
    lme1 = fitlme(alldata_table_SameHouse_roi,formula);
    formula = 'z~ entryRoom + (1|subj)';
    lme2 = fitlme(alldata_table_SameHouse_roi,formula);
    results = compare(lme1,lme2);
    disp('SH room') % 
    disp(results.pValue)
    mixed_effect_model_SH_p(roi_i) = results.pValue(2);
    figure(f0); axes('position',pos(roi_ind,:))
    pval_room(roi_i,:) = plot_bar_panel3(Z,xlabel_str, barcolor);
    axis([0,4,0,YTICKS(2)])
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)
 
    POSITION = get(gca,'position');
    
    
    XLIM = get(gca,'xlim');
    YLIM = get(gca,'Ylim');
    text(XLIM(1),YLIM(2)+0.5*diff(YLIM),{ sprintf('%s',regexprep(roi_str{roi_i},'_',' ')),...
                                             sprintf('SH mixed-effect: P = %1.2e',mixed_effect_model_SH_p(roi_i))} );
end


%%
title_str1 = 'fig2b_diffHouse';
f0 = newA4figure(title_str1);
pval_room = [];
roi_list = [1,2];
xlabel_str = {'entryRoom','innerRoom','cross room'};


for roi_ind = 1:length(roi_list)
    
    roi_i = roi_list(roi_ind);
    
    % restrict to PS of analyzed ROI - 
    row_num = find(ismember(alldata_table_house.roi,roi_str{roi_i}));
    alldata_table_house_roi = alldata_table_house(row_num,:);

    % restrict for trials from different videos exploring the same house - 
    row_num = find(ismember(alldata_table_house_roi.condition,condition_val{1}));
    alldata_table_DiffHouse_roi = alldata_table_house_roi(row_num,:);
    
    if roi_ind == 1
        YTICKS = [0 4]*10^-3;
    elseif roi_ind == 2
        YTICKS = [0 8]*10^-3;
    elseif roi_ind == 3
        YTICKS = [0 ,3.5]*10^-2;
    elseif roi_ind == 4
        YTICKS = [0 ,10]*10^-3;
    end
    
     
    clear Z
    rows = find( ismember(alldata_table_DiffHouse_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_DiffHouse_roi.z(rows);
    rows = find( ismember(alldata_table_DiffHouse_roi.entryRoom,'innerRoom'));
    Z{1,2}=  alldata_table_DiffHouse_roi.z(rows);
    rows = find( ismember(alldata_table_DiffHouse_roi.entryRoom,'otherLoc'));
    Z{1,3}=  alldata_table_DiffHouse_roi.z(rows);
   
    y = [Z{1,1}',Z{1,2}',Z{1,3}'];
    g = [1*ones(1,length(Z{1,1})),2*ones(1,length(Z{1,2})),3*ones(1,length(Z{1,3}))];
    [p_anovan(2+roi_i),tbl,stats] = anovan(y,{g});
    figure; multcompare(stats)
    
    formula = 'z~ 1 +(1|subj)';
    lme1 = fitlme(alldata_table_DiffHouse_roi,formula);
    formula = 'z~ entryRoom + (1|subj)';
    lme2 = fitlme(alldata_table_DiffHouse_roi,formula);
    results = compare(lme1,lme2);
    disp('SH room') % 
    disp(results.pValue)
    mixed_effect_model_SH_p(roi_i) = results.pValue(2);
    figure(f0); axes('position',pos(roi_ind,:))
    pval_room(roi_i,:) = plot_bar_panel3(Z,xlabel_str, barcolor);
    axis([0,4,0,YTICKS(2)])
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)
 
    POSITION = get(gca,'position');
    
    
    XLIM = get(gca,'xlim');
    YLIM = get(gca,'Ylim');
    text(XLIM(1),YLIM(2)+0.5*diff(YLIM),{ sprintf('%s',regexprep(roi_str{roi_i},'_',' ')),...
                                             sprintf('SH mixed-effect: P = %1.2e',mixed_effect_model_SH_p(roi_i))} );
end


%% STATS when aggregating both houses together
% sub panel - all houses - entryRoom
title_str1 = 'fig2b_supData1';
f0 = newA4figure(title_str1);

roi_list = [1,2];
for roi_ind = 1:length(roi_list)
    
    roi_i = roi_list(roi_ind);
    
    % restrict to PS of analyzed ROI - 
    row_num = find(ismember(alldata_table_house.roi,roi_str{roi_i}));
    alldata_table_house_roi = alldata_table_house(row_num,:);

    if roi_ind == 1
        YTICKS = [0 4]*10^-3;
    elseif roi_ind == 2
        YTICKS = [0 8]*10^-3;
    elseif roi_ind == 3
        YTICKS = [0 ,3.5]*10^-2;
    elseif roi_ind == 4
        YTICKS = [0 ,10]*10^-3;
    end
    
    clear Z
    rows = find( ismember(alldata_table_house_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_house_roi.z(rows);
    rows = find( ismember(alldata_table_house_roi.entryRoom,'innerRoom'));
    Z{1,2}=  alldata_table_house_roi.z(rows);
    
    axes('position',pos(roi_ind,:))
    pval_room(roi_i) = plot_bar_panel(Z,xlabel_str, barcolor);
    set(gca,'xticklabels',[])
    
    formula = 'z~ 1 +(1|subj)';
    lme1 = fitlme(alldata_table_house_roi,formula);
    formula = 'z~ entryRoom + (1|subj)';
    lme2 = fitlme(alldata_table_house_roi,formula);
    results = compare(lme1,lme2);
    disp('room') % 
    disp(results.pValue)
    mixed_effect_model_p(roi_i) = results.pValue(2);
    
    
    % For entry-room items - is house info reflected?
    entryRoom_val = unique(alldata_table_house_roi.entryRoom);
    row_num = find(ismember(alldata_table_house_roi.entryRoom,entryRoom_val{1}));
    alldata_table_entryRoom = alldata_table_house_roi(row_num,:);
    unique(alldata_table_entryRoom.condition)
    
    formula = 'z~ 1 +(1|subj)';
    lme1 = fitlme(alldata_table_entryRoom,formula);
    formula = 'z~ condition + (1|subj)';
    lme2 = fitlme(alldata_table_entryRoom,formula);
    results = compare(lme1,lme2);
    disp('All item-pairs from entry-room - testing for house-identity')
    disp(results.pValue) % P


    %title(sprintf('%s',regexprep(roi_str{roi_i},'_',' ')))
    axis([0,3,0,inf])
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)

    XLIM = get(gca,'xlim');
    YLIM = get(gca,'Ylim');
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)
    text(XLIM(1),YLIM(2)+0.5*diff(YLIM),{ sprintf('%s',regexprep(roi_str{roi_i},'_',' ')),...
                                             sprintf('SH/DH mixed: P = %1.2e',mixed_effect_model_p(roi_i)) } );
end


%% STATS when aggregating both houses together
% sub panel - all houses - entryRoom
title_str1 = 'fig2b_supData2';
f0 = newA4figure(title_str1);

roi_list = [1,2];
for roi_ind = 1:length(roi_list)
    
    roi_i = roi_list(roi_ind);
    
    % restrict to PS of analyzed ROI - 
    row_num = find(ismember(alldata_table_house.roi,roi_str{roi_i}));
    alldata_table_house_roi = alldata_table_house(row_num,:);

    % restrict for trials from different videos exploring *different* houses - 
    % Motivation - 
    row_num = find(ismember(alldata_table_house_roi.condition,condition_val{1}));
    alldata_table_DiffHouse_roi = alldata_table_house_roi(row_num,:);
  
    
    if roi_ind == 1
        YTICKS = [0 4]*10^-3;
    elseif roi_ind == 2
        YTICKS = [0 8]*10^-3;
    elseif roi_ind == 3
        YTICKS = [0 ,3.5]*10^-2;
    elseif roi_ind == 4
        YTICKS = [0 ,10]*10^-3;
    end
    
    clear Z
    rows = find( ismember(alldata_table_DiffHouse_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_DiffHouse_roi.z(rows);
    rows = find( ~ismember(alldata_table_DiffHouse_roi.entryRoom,'entryRoom'));
    Z{1,2}=  alldata_table_DiffHouse_roi.z(rows);
    
    axes('position',pos(roi_ind,:))
    pval_room(roi_i) = plot_bar_panel(Z,xlabel_str, barcolor);
    set(gca,'xticklabels',[])
    
    formula = 'z~ 1 +(1|subj)';
    lme1 = fitlme(alldata_table_DiffHouse_roi,formula);
    formula = 'z~ entryRoom + (1|subj)';
    lme2 = fitlme(alldata_table_DiffHouse_roi,formula);
    results = compare(lme1,lme2);
    disp('room') % 
    disp(results.pValue)
    mixed_effect_model_p(roi_i) = results.pValue(2);
    
    %title(sprintf('%s',regexprep(roi_str{roi_i},'_',' ')))
    axis([0,3,0,inf])
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)

    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)
    text(XLIM(1),YLIM(2)+0.5*diff(YLIM),{ sprintf('%s',regexprep(roi_str{roi_i},'_',' ')),...
                                             sprintf('DH mixed: P = %1.2e',mixed_effect_model_p(roi_i)) } );
end

%%
PrintActiveFigs('E:\Dropbox\RanganathLab\spatcon_2\figures\LINKS')

%% Alternative version
x_w = 0.08;
pos(1,:) = [0.1, 0.1, x_w, 0.15];
pos(2,:) = [0.27, 0.1, x_w, 0.15];
pos(3,:) = [0.44, 0.1, x_w, 0.15];
pos(4,:) = [0.7, 0.1, x_w, 0.15];

title_str1 = 'fig3b';
f0 = newA4figure(title_str1);
pval_room = [];
roi_list = [1,2];
for roi_ind = 1:length(roi_list)
    
    roi_i = roi_list(roi_ind);
    
    row_num = find(ismember(alldata_table_house.roi,roi_str{roi_i}));
    alldata_table_house_roi = alldata_table_house(row_num,:);

    row_num = find(ismember(alldata_table_house_roi.condition,condition_val{2}));
    alldata_table_SameHouse_roi = alldata_table_house_roi(row_num,:);
    
    if roi_ind == 1
        YTICKS = [0 4]*10^-3;
    elseif roi_ind == 2
        YTICKS = [0 8]*10^-3;
    elseif roi_ind == 3
        YTICKS = [0 ,3.5]*10^-2;
    elseif roi_ind == 4
        YTICKS = [0 ,10]*10^-3;
    end
    cmap = brewermap(5,'set1');
    barcolor(1,:) = cmap(5,:);
    barcolor(2,:) = cmap(3,:);
    
    % main panel - same house entryRoom
    clear Z
    rows = find( ismember(alldata_table_SameHouse_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_SameHouse_roi.z(rows);
    rows = find( ~ismember(alldata_table_SameHouse_roi.entryRoom,'entryRoom'));
    Z{1,2}=  alldata_table_SameHouse_roi.z(rows);
    formula = 'z~ 1 +(1|subj)';
    lme1 = fitlme(alldata_table_SameHouse_roi,formula);
    formula = 'z~ entryRoom + (1|subj)';
    lme2 = fitlme(alldata_table_SameHouse_roi,formula);
    results = compare(lme1,lme2);
    disp('SH room') % 
    disp(results.pValue)
    mixed_effect_model_SH_p(roi_i) = results.pValue(2);
    figure(f0); axes('position',pos(roi_ind,:))
    xlabel_str = {'entryRoom','other'};
    pval_room(roi_i) = plot_bar_panel(Z,xlabel_str,barcolor);
    axis([0,3,0,inf])
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)
 
    POSITION = get(gca,'position');
    
    % sub panel - all houses - entryRoom
    clear Z
    rows = find( ismember(alldata_table_house_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_house_roi.z(rows);
    rows = find( ~ismember(alldata_table_house_roi.entryRoom,'entryRoom'));
    Z{1,2}=  alldata_table_house_roi.z(rows);
    
    axes('position',[POSITION(1)+.85*POSITION(3), POSITION(2)+.87*POSITION(4),POSITION(3)*0.5,POSITION(4)*0.5])
    pval_room(roi_i) = plot_bar_panel(Z,xlabel_str, barcolor);
    set(gca,'xticklabels',[])
    
    formula = 'z~ 1 +(1|subj)';
    lme1 = fitlme(alldata_table_house_roi,formula);
    formula = 'z~ entryRoom + (1|subj)';
    lme2 = fitlme(alldata_table_house_roi,formula);
    results = compare(lme1,lme2);
    disp('room') % 
    disp(results.pValue)
    mixed_effect_model_p(roi_i) = results.pValue(2);
    
    %title(sprintf('%s',regexprep(roi_str{roi_i},'_',' ')))
    axis([0,3,0,inf])
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)

   
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)
    
    XLIM = get(gca,'xlim');
    YLIM = get(gca,'Ylim');
     text(XLIM(1),YLIM(2)+0.5*diff(YLIM),{ sprintf('%s',regexprep(roi_str{roi_i},'_',' ')),...
                                             sprintf('SH mixed: P = %1.2e',mixed_effect_model_SH_p(roi_i)),...
                                             sprintf('SH/DH mixed: P = %1.2e',mixed_effect_model_p(roi_i)) } );
end

PrintActiveFigs('E:\Dropbox\RanganathLab\spatcon_2\figures\LINKS')



%%
title_str1 = 'fig3c';
f0 = newA4figure(title_str1);
pval_room = [];
roi_list = [1 2];
for roi_ind = 1:length(roi_list)
    
    roi_i = roi_list(roi_ind);
    
    clear alldata_table_house_roi
    row_num = find(ismember(alldata_table_house.roi,roi_str{roi_i}));
    alldata_table_house_roi = alldata_table_house(row_num,:);
            
    rows1 = find(ismember(alldata_table_house_roi.entryRoom,'entryRoom'));
    alldata_table_house_roi_entryRoom = alldata_table_house_roi(rows1,:);
        
    formula = 'z~ 1 +(1|subj)';
    lme1 = fitlme(alldata_table_house_roi_entryRoom,formula);
    formula = 'z~ condition + (1|subj)';
    lme2 = fitlme(alldata_table_house_roi_entryRoom,formula);
    results = compare(lme1,lme2);
    disp('same/diff house') % 
    disp(results.pValue)
    mixed_effect_model_entryRoom_p(roi_i) = results.pValue(2);

    clear Z
    rows = find(ismember(alldata_table_house_roi_entryRoom.condition,'differentVideo_sameHouse'));
    Z{1,1} =  alldata_table_house_roi_entryRoom.z(rows);
    
    rows = find(ismember(alldata_table_house_roi_entryRoom.condition,'differentHouse_eitherRoom'));
    Z{1,2} =  alldata_table_house_roi_entryRoom.z(rows);
    
    [pval_test,h] = ranksum(Z{1,1},Z{1,2})

    figure(f0); axes('position',pos(roi_ind,:))
    xlabel_str = {'same house','diff house'};
    pval_room(roi_i) = plot_bar_panel(Z,xlabel_str);
 
    if roi_ind == 1
        YTICKS = [0 4]*10^-3;
    elseif roi_ind == 2
        YTICKS = [0 7]*10^-3;
    elseif roi_ind == 3
        YTICKS = [0 ,15]*10^-3;
    end
    set(gca,'ylim',YTICKS,'ytick',YTICKS,'yticklabel',YTICKS)
    
    XLIM = get(gca,'xlim');
    YLIM = get(gca,'Ylim');
    text(XLIM(1),YLIM(2)+0.5*diff(YLIM),{sprintf('%s',regexprep(roi_str{roi_i},'_',' ')),...
                                            sprintf('bootstrapping: P = %2.4f',p_space_BS(roi_i)),...
                                            sprintf('Wilcoxson: P = %2.3f',pval_room(roi_i))});

end
PrintActiveFigs('E:\Dropbox\RanganathLab\spatcon_2\figures\LINKS')



%% RUN BOOTSRAPPING TO VERIFY RSA BAR COMPARISON
pos(1,:) = [0.1, 0.1, 0.12, 0.15];
pos(2,:) = [0.25, 0.1, 0.12, 0.15];
pos(3,:) = [0.4, 0.1, 0.12, 0.15];
pos(4,:) = [0.1, 0.4, 0.12, 0.15];
pos(5,:) = [0.25, 0.4, 0.12, 0.15];
pos(6,:) = [0.4, 0.4, 0.12, 0.15];

title_str1 = 'room_code_ns';
title_str2 = 'house_and_room_code_ns';

f0 = newA4figure(title_str1);
f1 = newA4figure(title_str2);

outputFolder = 'E:\spatcon_hdz\data_p\corrAnalysis_permutation\entryRoom_fig2B_bars';
roi_list = [1,2,4,5:7];
run_bootstrapping = 0;
for roi_ind = 1:length(roi_list);
    
    roi_i = roi_list(roi_ind);
    
    row_num = find(ismember(alldata_table_house.roi,roi_str{roi_i}));
    alldata_table_house_hem_roi = alldata_table_house(row_num,:);
    
    % entryRoom
    clear Z
    rows = find( ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_house_hem_roi.z(rows);
    rows = find( ~ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
    Z{1,2}=  alldata_table_house_hem_roi.z(rows);
    
    if run_bootstrapping
        N_BS_trials = 10000;
        for ii_rb = 1:N_BS_trials
            x = Z{1,1};
            k = length(x);
            y = datasample(Z{1,2},k);
            m(ii_rb) = median(x)-median(y);
        end
        p_BS = sum(m > 0)/length(m);
        save(fullfile(outputFolder,sprintf('allZ_struct_roi%d',roi_i)),'Z','p_BS','N_BS_trials','m')
    end
    
    figure(f0); axes('position',pos(roi_ind,:))
    xlabel_str = {'entryRoom','other'};
    plot_2by2_distributions(Z,xlabel_str)
    [pval_room(hem_i, roi_i),h] = ranksum(Z{1,1},Z{1,2});
    
    title(sprintf('%s',regexprep(roi_str{roi_i},'_',' ')))
    axis([0,3,0,inf])
    
   % Limiting to same house - comparing entryRoom to others   
    row_num = find(~ismember(alldata_table_house_hem_roi.condition,'differentVideo_sameHouse'));
    alldata_table_same_house_hem_roi = alldata_table_house_hem_roi(row_num,:);

    clear Z
    rows = find( ismember(alldata_table_same_house_hem_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_same_house_hem_roi.z(rows);
    rows = find( ~ismember(alldata_table_same_house_hem_roi.entryRoom,'entryRoom'));
    Z{1,2}=  alldata_table_same_house_hem_roi.z(rows);
    
    if run_bootstrapping
        N_BS_trials = 10000;
        for ii_rb = 1:N_BS_trials
            x = Z{1,1};
            k = length(x);
            y = datasample(Z{1,2},k);
            m(ii_rb) = median(x)-median(y);
        end
        p_SH_BS = sum(m > 0)/length(m);
        save(fullfile(outputFolder,sprintf('BS_SH_allZ_struct_roi%d',roi_i)),'Z','p_SH_BS','N_BS_trials','m')
    end
    
    figure(f0); axes('position',pos(roi_ind,:))
    xlabel_str = {'SH entryRoom','other'};
    plot_2by2_distributions(Z,xlabel_str)
    set(gca, 'XTickLabelRotation',30)
    [pval_SH_entryRoom(roi_i),h] = ranksum(Z{1,1},Z{1,2});
    
    title(sprintf('%s',regexprep(roi_str{roi_i},'_',' ')))
    axis([0,3,0,inf])
    
    
    % Figure 2
    % same house same room
    clear Z
    rows = find(ismember(alldata_table_house_hem_roi.condition,'differentVideo_sameHouse') &...
        ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
    Z{1,1} =  alldata_table_house_hem_roi.z(rows);
    
    rows = find(~ismember(alldata_table_house_hem_roi.condition,'differentHouse_eitherRoom') |...
        ismember(alldata_table_house_hem_roi.entryRoom,'entryRoom'));
    Z{1,2} =  alldata_table_house_hem_roi.z(rows);
    
    outputFolder = 'E:\spatcon_hdz\data_p\corrAnalysis_permutation\entryRoom_SH_DH_fig2C';
    if run_bootstrapping
        N_BS_trials = 10000;
        for ii_rb = 1:N_BS_trials
            x = Z{1,1};
            k = length(x);
            y = datasample(Z{1,2},k);
            m(ii_rb) = median(x)-median(y);
        end
        p_BS = sum(m > 0)/length(m);
        save(fullfile(outputFolder,sprintf('space_allZ_struct_roi%d',roi_i)),'Z','p_BS','N_BS_trials','m')
    end
    
    figure(f1); axes('position',pos(roi_ind,:))
    xlabel_str = {'SH/R1','DH/R1'};
    plot_2by2_distributions(Z,xlabel_str)
    title(sprintf('%s',regexprep(roi_str{roi_i},'_',' ')))

    [pval_space(hem_i, roi_i),h] = ranksum(Z{1,1},Z{1,2});
    axis([0,3,0,inf])
        
end

PrintActiveFigs('E:\Dropbox\RanganathLab\spatcon_2\figures\LINKS');
%%

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
        formula = 'z~ spaceCond + (1|subj)';
        lme2 = fitlme(alldata_table_house_hem_roi,formula);
        results = compare(lme1,lme2);
        test_str{test_i} = 'spaceCond - DV - SH/SR vs others';
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

        % Test items in entry room vs all others - for Fig 3B
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

    end
end

if run_bootstrapping
    % Summarize and plot the results
    outputDir = 'E:\spatcon_hdz\data_p\corrAnalysis_permutation\entryRoom_fig2B_bars';
    a = dir(fullfile(outputDir,'allZ_struct_roi*'));
    p_BS= NaN(1,7);
    for ii = 1:length(a)
        roi_num = str2num(a(ii).name(end-4));
        mm = matfile(fullfile(outputDir,a(ii).name));
        p_BS(roi_num) = mm.p_BS;
        if  p_BS(ii) > 0.2;
            p_BS(roi_num) = 1 - mm.p_BS;
        end
    end
    
    outputDir = 'E:\spatcon_hdz\data_p\corrAnalysis_permutation\entryRoom_SH_DH_fig2C';
    a = dir(fullfile(outputDir,'space_allZ_struct_roi*'));
    p_space_BS= NaN(1,7);
    for ii = 1:length(a)
        roi_num = str2num(a(ii).name(end-4));
        mm = matfile(fullfile(outputDir,a(ii).name));
        p_space_BS(roi_num) = mm.p_BS;
        if  p_space_BS(roi_num) > 0.2;
            p_space_BS(roi_num) = 1 - mm.p_BS;
        end
    end
    
end

end % func


