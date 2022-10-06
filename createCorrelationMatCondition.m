% This code uses data-files created by corrExtractionALL_publishedData.m
% This function creates hem/roi matrices specific to every condition
% sameHouse_sameVideo
% sameHouse_diffVideo
% sameHouseRoom_diffVideo
% sameHouseRoomLoc_diffVideo

% Before running this - also complete - 
% stimuli extraction per subject - stimuliExtraction.m and loadEncodingTesting.m

function createCorrelationMatCondition

% List order is determined by corrExtractionALL_publishedData
% rois_of_most_interest= {'CA1_body', 'CA2_3_DG_body',...
%    'ERC', 'whole_hippo', 'subiculum','35_36'};

% subject sepcific onset files are extracted using stimuliExtraction.m
onsetFiles_root = 'E:\spatcon_hdz\data_p\onsetFiles\';

for roiSet_i = 1:2
    
    if roiSet_i == 1
        mm = matfile('E:\spatcon_hdz\data_p\corrAnalysis\pattern_corr_zscored_all_subjects.mat'); % extracted by corrExtractionALL_localData/corrExtractionALL_publishedData
        object_info_all_rois = mm.object_info_all_rois;
        pattern_corr_all_rois = mm.pattern_corr_all_rois;
        savefile_Str = '';
    elseif roiSet_i == 2
        mm = matfile('E:\spatcon_hdz\data_p\corrAnalysis\pattern_corr_manualTracing_zscored_all_subjects.mat'); % extracted by corrExtractionALL_localData/corrExtractionALL_manuallyTracedData
        object_info_all_rois = mm.object_info_all_rois;
        pattern_corr_all_rois = mm.pattern_corr_all_rois;
        savefile_Str = 'manualTraces_';
    end
    [N,H,R] = size(pattern_corr_all_rois);


for cond_i = [4] 
    
    COLLECT_ALL = 0;
    COLLECT_Rhit = 0;
    COLLECT_Rhit_HC = 0;
    COLLECT_Rhit_HC_RC = 0;
    COLLECT_errors = 0;
    COLLECT_roomLoc = 0;
    COLLECT_exactRoomLoc = 0;
    COLLECT_Rhit_HC_roomLoc = 0;
    COLLECT_post_RoomErrors = 0;
    COLLECT_mixed_encoding_status = 0;
    COLLECT_Rhit_HC_encRoomLoc = 0;
    
    disp(['cond ',num2str(cond_i)])
    switch cond_i
        case 1
            COLLECT_ALL = 1;
        case 2
            COLLECT_Rhit = 1;
        case 3
            COLLECT_Rhit_HC = 1;
        case 4
            COLLECT_Rhit_HC_RC = 1;
        case 5
            COLLECT_errors = 1;
        case 6
            COLLECT_roomLoc = 1; % we want cases where Rhit = 1, encoding-RC = 1
        case 7
            COLLECT_exactRoomLoc = 1;
        case 8
            COLLECT_Rhit_HC_roomLoc = 1; 
        case 9
            COLLECT_post_RoomErrors = 1;
        case 10
            COLLECT_mixed_encoding_status = 1;
        case 11
            COLLECT_Rhit_HC_roomLocErrors = 1;
        end
    
    for hem_i = 1:H
        
        if hem_i == 1
            hemi_label = 'ashs_left';
        else
            hemi_label = 'ashs_right';
        end
        
        for roi_i = 1:R
            
            % condition - same house same video
            sameHouse_sameVideo = createStruct();
            % condition - same house diff video
            sameHouse_diffVideo  = createStruct();
            % condition - same house/room diff video
            sameHouseRoom_diffVideo  = createStruct();
            % condition - same house/room/loc in room diff video
            sameHouseRoomLoc_diffVideo  = createStruct();
            % condition - room #1 in diff video
            sameRoomLoc  = createStruct();
            % condition - room #1 in diff video
            diffRoomLoc  = createStruct();
            
            
            % Controls -
            % condition - diff house, diff video
            diffHouse_diffVideo  = createStruct();
            % condition - all error trials
            allTrials_Rmiss  = createStruct();
            
            
            for subj = 1:N
                
                if roiSet_i == 2 && subj == 28; continue; end %%%%% DEBUG - matrices for s028 are in the wrong size !!!!

                if isempty(object_info_all_rois{subj}); continue; end
                if isempty(pattern_corr_all_rois{subj,hem_i, roi_i}); continue; end
                
                
                folderName = fullfile(onsetFiles_root, sprintf('s%03d',subj));
                try
                    stimuliFile = fullfile(folderName,sprintf('stimuliRespList_s%03d.mat',subj));
                    mm = matfile(stimuliFile);
                    objectID = mm.objectID;
                    REF_ID = mm.REF_ID;
                    
                catch
                    error('onset file not found')
                end
                
                
                
                for ii_o = 1:length(object_info_all_rois{subj}.object_num)
                    for jj_o = 1:length(object_info_all_rois{subj}.object_num)
                        
                        if ii_o == jj_o; continue; end
                        if isnan(object_info_all_rois{subj}.object_num(ii_o)); continue; end
                        if isnan(object_info_all_rois{subj}.object_num(jj_o)); continue; end
                        
                        id1 = find(objectID(:,REF_ID.object_num) == object_info_all_rois{subj}.object_num(ii_o));
                        ID1 = objectID(id1,:);
                        id2 = find(objectID(:,REF_ID.object_num) == object_info_all_rois{subj}.object_num(jj_o));
                        ID2 = objectID(id2,:);
                        
                        if (isnan(pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o))); continue; end
                        if (pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o) == 0); continue; end
                        
                        % same block correlations are removed - high correlations regardless
                        % of item
                        if ID1(REF_ID.BlockID) == ID2(REF_ID.BlockID); continue; end
                                                                        
                        if ~COLLECT_ALL && COLLECT_Rhit
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                        elseif ~COLLECT_ALL && COLLECT_Rhit_HC
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if ID1(REF_ID.HC) == 0 || ID2(REF_ID.HC) == 0; continue; end % remove items that weren't remembered in the correct house
                        elseif ~COLLECT_ALL && COLLECT_Rhit_HC_RC
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if ID1(REF_ID.RC) == 0 || ID2(REF_ID.RC) == 0; continue; end % remove items that weren't remembered in the correct house
                        elseif COLLECT_roomLoc % Rhit = 1, encoding-RC = 1 
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if ID1(REF_ID.LocInRoomRespBallPark) == 0 || ID2(REF_ID.LocInRoomRespBallPark) == 0; continue; end % remove items that weren't remembered in the correct house
                        elseif COLLECT_errors  % Rhit = 1, encoding-RC = 0
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if ID1(REF_ID.LocInRoomRespBallPark) == 1 || ID2(REF_ID.LocInRoomRespBallPark) == 1; continue; end % remove items that were remembered well
                        elseif COLLECT_mixed_encoding_status  % Rhit = 1, encoding-RC = 0
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if ID1(REF_ID.LocInRoomRespBallPark) == ID2(REF_ID.LocInRoomRespBallPark); continue; end % remove couples with identical status (both encoded/both not encoded)
                        elseif COLLECT_exactRoomLoc
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if (ID1(REF_ID.LocInRoomResp) ~= ID1(REF_ID.locationInRoom)) || (ID2(REF_ID.LocInRoomResp) ~= ID2(REF_ID.locationInRoom)) ; continue; end
                        elseif COLLECT_Rhit_HC_roomLoc % we want cases where Rhit = 1, HC = 1 , encoding-RC = 1
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if ID1(REF_ID.HC) == 0 || ID2(REF_ID.HC) == 0; continue; end % remove items that weren't remembered in the correct house
                            if ID1(REF_ID.LocInRoomRespBallPark) == 0 || ID2(REF_ID.LocInRoomRespBallPark) == 0; continue; end % remove items that weren't encoded in the correct room
                        elseif COLLECT_post_RoomErrors % we want cases where Rhit = 1, HC = 1 but RC = 0
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if ID1(REF_ID.HC) == 0 || ID2(REF_ID.HC) == 0; continue; end % remove items that weren't remembered in the correct house
                            if ID1(REF_ID.RC) ~= 0 || ID2(REF_ID.RC) ~= 0; continue; end
                        elseif COLLECT_Rhit_HC_roomLocErrors % we want cases where Rhit = 1, HC = 1 , encoding-RC = 0
                            if ID1(REF_ID.Rhit) == 0 || ID2(REF_ID.Rhit) == 0; continue; end % remove items that weren't remembered
                            if ID1(REF_ID.HC) == 0 || ID2(REF_ID.HC) == 0; continue; end % remove items that weren't remembered in the correct house
                            if ID1(REF_ID.LocInRoomRespBallPark) == 1 || ID2(REF_ID.LocInRoomRespBallPark) == 1; continue; end % remove items that were remembered well
                        end
                        
                        sameHouse = 0; sameRoom = 0; sameVideo = 0; sameLoc = 0; boundaryLoc = 0; entryRoom = 0;
                        if ID1(REF_ID.house) == ID2(REF_ID.house); sameHouse = 1; end
                        if ID1(REF_ID.room) == ID2(REF_ID.room); sameRoom = 1; end
                        if ( ID1(REF_ID.room) == 1 && ID2(REF_ID.room) == 1 ); entryRoom = 1; end
                        if ID1(REF_ID.video) == ID2(REF_ID.video); sameVideo = 1; end
                        
                        
                        
                        % locations in room are 1-8
                        if ( sum(ismember([ID1(REF_ID.locationInRoom),ID2(REF_ID.locationInRoom)],[1,2]))== 2 ||...
                             sum(ismember([ID1(REF_ID.locationInRoom),ID2(REF_ID.locationInRoom)],[2,3]))== 2 ||...     
                             sum(ismember([ID1(REF_ID.locationInRoom),ID2(REF_ID.locationInRoom)],[3,4]))== 2 ||...
                             sum(ismember([ID1(REF_ID.locationInRoom),ID2(REF_ID.locationInRoom)],[5,6]))== 2 ||...
                             sum(ismember([ID1(REF_ID.locationInRoom),ID2(REF_ID.locationInRoom)],[6,7]))== 2 ||...
                             sum(ismember([ID1(REF_ID.locationInRoom),ID2(REF_ID.locationInRoom)],[7,8]))== 2 )
                             sameLoc = 1; 
                        end
                        
                        % near doors 
                        if ( sum(ismember([ID1(REF_ID.locationInRoom),ID2(REF_ID.locationInRoom)],[1,2,5,6]))==2 )
                             boundaryLoc = 1; 
                        end
                        
                        % same room, regardless of house/vid info
                        if entryRoom && boundaryLoc && ~sameVideo
                            sameRoomLoc = updateStruct(sameRoomLoc, subj,...
                                object_info_all_rois{subj}.pattern_ids_all_runs{ii_o},...
                                object_info_all_rois{subj}.pattern_ids_all_runs{jj_o},...
                                object_info_all_rois{subj}.object_num(ii_o),...
                                object_info_all_rois{subj}.object_num(jj_o),...
                                pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o));
                        elseif ~sameRoom && ~sameVideo
                            diffRoomLoc = updateStruct(diffRoomLoc, subj,...
                                object_info_all_rois{subj}.pattern_ids_all_runs{ii_o},...
                                object_info_all_rois{subj}.pattern_ids_all_runs{jj_o},...
                                object_info_all_rois{subj}.object_num(ii_o),...
                                object_info_all_rois{subj}.object_num(jj_o),...
                                pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o));
                        end

                        if sameHouse && sameVideo
                            sameHouse_sameVideo = updateStruct(sameHouse_sameVideo, subj,...
                                object_info_all_rois{subj}.pattern_ids_all_runs{ii_o},...
                                object_info_all_rois{subj}.pattern_ids_all_runs{jj_o},...
                                object_info_all_rois{subj}.object_num(ii_o),...
                                object_info_all_rois{subj}.object_num(jj_o),...
                                pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o));
                        end
                        
                        if sameHouse && ~sameVideo
                            sameHouse_diffVideo = updateStruct(sameHouse_diffVideo, subj,...
                                object_info_all_rois{subj}.pattern_ids_all_runs{ii_o},...
                                object_info_all_rois{subj}.pattern_ids_all_runs{jj_o},...
                                object_info_all_rois{subj}.object_num(ii_o),...
                                object_info_all_rois{subj}.object_num(jj_o),...
                                pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o));
                        end
                        
                        if ~sameHouse && ~sameVideo
                            diffHouse_diffVideo = updateStruct(diffHouse_diffVideo, subj,...
                                object_info_all_rois{subj}.pattern_ids_all_runs{ii_o},...
                                object_info_all_rois{subj}.pattern_ids_all_runs{jj_o},...
                                object_info_all_rois{subj}.object_num(ii_o),...
                                object_info_all_rois{subj}.object_num(jj_o),...
                                pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o));
                        end
                        
                        if sameHouse && sameRoom && ~sameVideo
                            sameHouseRoom_diffVideo = updateStruct(sameHouseRoom_diffVideo, subj,...
                                object_info_all_rois{subj}.pattern_ids_all_runs{ii_o},...
                                object_info_all_rois{subj}.pattern_ids_all_runs{jj_o},...
                                object_info_all_rois{subj}.object_num(ii_o),...
                                object_info_all_rois{subj}.object_num(jj_o),...
                                pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o));
                        end
                        
                        
                        if sameHouse && entryRoom && boundaryLoc && ~sameVideo
                            sameHouseRoomLoc_diffVideo = updateStruct(sameHouseRoomLoc_diffVideo, subj,...
                                object_info_all_rois{subj}.pattern_ids_all_runs{ii_o},...
                                object_info_all_rois{subj}.pattern_ids_all_runs{jj_o},...
                                object_info_all_rois{subj}.object_num(ii_o),...
                                object_info_all_rois{subj}.object_num(jj_o),...
                                pattern_corr_all_rois{subj,hem_i, roi_i}(ii_o,jj_o));
                        end
                        
                    end
                end
            end
            
            % save structs for hem/roi
            sameHouse_sameVideo_all = sameHouse_sameVideo;
            sameHouse_diffVideo_all = sameHouse_diffVideo;
            sameHouseRoom_diffVideo_all = sameHouseRoom_diffVideo;
            sameHouseRoomLoc_diffVideo_all  = sameHouseRoomLoc_diffVideo;
            diffHouse_diffVideo_all = diffHouse_diffVideo;
            sameRoomLoc_all = sameRoomLoc;
            diffRoomLoc_all = diffRoomLoc;
            
            if COLLECT_ALL
                filename = sprintf('pattern_corr_zscored_conditions_%shem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_Rhit
                filename = sprintf('pattern_corr_zscored_conditions_%sRhit_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_Rhit_HC
                filename = sprintf('pattern_corr_zscored_conditions_%sRhit_HC_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_Rhit_HC_RC
                filename = sprintf('pattern_corr_zscored_conditions_%sRhit_HC_RC_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_errors
                filename = sprintf('pattern_corr_zscored_conditions_%serrors_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_roomLoc
                filename = sprintf('pattern_corr_zscored_conditions_%sRlocC_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_exactRoomLoc
                 filename = sprintf('pattern_corr_zscored_conditions_%sRlocC_exact_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_post_RoomErrors
                 filename = sprintf('pattern_corr_zscored_conditions_%spost_RoomErrors_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_Rhit_HC_roomLoc
                 filename = sprintf('pattern_corr_zscored_conditions_%sRhit_HC_RlocC_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_mixed_encoding_status
                filename = sprintf('pattern_corr_zscored_conditions_%sRhit_mixed_RlocC_error_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);
            elseif COLLECT_Rhit_HC_roomLocErrors
                 filename = sprintf('pattern_corr_zscored_conditions_%sRhit_HC_RlocError_hem%d_roi_%d.mat',savefile_Str,hem_i,roi_i);

            end
            
            disp(['saving ',filename])

            save(fullfile('E:\spatcon_hdz\data_p\corrAnalysis\',filename),'sameHouse_sameVideo_all',...
                'sameHouse_diffVideo_all','sameHouseRoom_diffVideo_all','sameHouseRoomLoc_diffVideo_all',...
                'diffHouse_diffVideo_all','allTrials_Rmiss','diffRoomLoc_all','sameRoomLoc_all','-v7.3')
            
            
            clear sameHouse_sameVideo_all sameHouse_diffVideo_all sameHouseRoom_diffVideo_all sameHouseRoomLoc_diffVideo_all diffHouse_diffVideo_all
            
        end
    end
    
end

end
end

function condStruct= createStruct()

condStruct.subj = [];
condStruct.col1 = {};
condStruct.col2 = {};
condStruct.objID1 = [];
condStruct.objID2 = [];
condStruct.fisherZ_r = [];

end

function condStruct = updateStruct(condStruct, subj, col1,col2, id1,id2, fz_r)

condStruct.subj(end+1) = subj;
condStruct.col1{end+1} = col1;
condStruct.col2{end+1} = col2;
condStruct.objID1(end+1) = id1;
condStruct.objID2(end+1) = id2;
condStruct.fisherZ_r(end+1) = fz_r;

end
