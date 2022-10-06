% Loading dat files of testing during encoding session, to grab the
% accuracy of spatial-location whithin a room, right after viewing the
% video

function loadEncodingTesting

outputFigures = 'E:\spatcon_hdz\data_p\figures\outputFolder';
PLOT = 0;

% subject sepcific onset files are extracted using stimuliExtraction.m
onsetFiles_root = 'E:\spatcon_hdz\data_p\onsetFiles\';
behaviorDataDir = 'E:\spatcon_hdz\data\raw_behavioral\';

subjects = 1:30;
exclude_subjects = [3:5,15:17,22]; % Excluded from original publication
subjects = subjects(~ismember(subjects,exclude_subjects));

allSubjectRating.Loc= NaN(1,30);
allSubjectRating.Rhit= NaN(1,30);
allSubjectRating.HC= NaN(1,30);
allSubjectRating.RC= NaN(1,30);

cnt = 0;
if PLOT
    f1 = newA4figure('correctRoomAssignment');
    f2 = newA4figure('correctSectorAssignment');
end

itemNumber = 1:240;
houseType = zeros(1,length(itemNumber));
videoNum= zeros(1,length(itemNumber));
locInRoom_GT= zeros(1,length(itemNumber));
room_GT= zeros(1,length(itemNumber));

subj_view = zeros(length(itemNumber),length(subjects));
subj_correct = nan(length(itemNumber),length(subjects));
subj_incorrect = nan(length(itemNumber),length(subjects));
LocInRoomResponse_all = zeros(length(itemNumber),length(subjects));

for subj_i = 1:length(subjects)
    subj = subjects(subj_i);
    
    
    folderName = fullfile(onsetFiles_root, sprintf('s%03d',subj));
    try
        stimuliFile = fullfile(folderName,sprintf('stimuliList_s%03d.mat',subj));
        mm = matfile(stimuliFile);
        objectID = mm.objectID;
        REF_ID = mm.REF_ID;
        
    catch
        warning('stimuliFile file not found')
        continue
    end
    
    
    behaviorDataDirSubj = fullfile(behaviorDataDir,sprintf('s%03d',subj));
    encodingResFile = fullfile(behaviorDataDirSubj, sprintf('ConABCD_objectEnc_s%03d.dat',subj));
    A = readtable(encodingResFile);
    
    correctLocation = A.CorrectLocation; % there are <2>X<8> options in this slot
    subjectResponseLocation = A.LocationResp; % there are 8 options (sectors) 
    objid_str = A.ObjectID;
    for ii_o = 1:length(objid_str)
        object_num(ii_o) = str2num(objid_str{ii_o}(end-2:end));
    end
    subj_view(object_num,subj_i) = true;
    
    REF_ID.LocInRoomResp = 13;
    REF_ID.LocInRoomRespBallPark = 14;
    REF_ID.LocInRoomRespExact = 15;

    LearnedSet = [];
    for ii = 1:length(objid_str)
        ind = find(object_num(ii) == objectID(:,REF_ID.object_num));
        LearnedSet(end+1) = ind;
        
        correctLocation_obj = correctLocation(ii);
        correctLocation_inRoom = mod(correctLocation(ii),10);
        LocInRoomResponse = str2num(subjectResponseLocation{ii}(1));
        
        if isempty(LocInRoomResponse) % no response
            objectID(ind,REF_ID.LocInRoomResp) = NaN;
            objectID(ind,REF_ID.LocInRoomRespBallPark) = NaN;
            continue
        end
        
        objectID(ind,REF_ID.LocInRoomResp) = LocInRoomResponse;
        
        room_id = objectID(ind,REF_ID.room);
        LocRoom_id = objectID(ind,REF_ID.locationInRoom);
        
        % sanity check
        if correctLocation ~= (room_id*10 + LocRoom_id)
            error('encoding file different from main onset file')
        end
        
        objectID(ind,REF_ID.LocInRoomResp) = NaN; % need to fix to account for 8 vs 4 locations in ground-truth vs response
        
        objectID(ind,REF_ID.LocInRoomRespBallPark) = (room_id == 1 && LocInRoomResponse <=4) || ...
                                                     (room_id == 2 && LocInRoomResponse > 4);
                                                 
        objectID(ind,REF_ID.LocInRoomRespExact) = NaN; % need to fix to account for 8 vs 4 locations in ground-truth vs response
        
        if (objectID(ind,REF_ID.LocInRoomRespBallPark))
            subj_correct(object_num(ii), subj_i) = true;
            subj_incorrect(object_num(ii), subj_i) = false;
        else
            subj_incorrect(object_num(ii), subj_i) = true;
            subj_correct(object_num(ii), subj_i) = false;
        end
        
        houseType(object_num(ii)) = objectID(ind,REF_ID.house);
        videoNum(object_num(ii)) = objectID(ind,REF_ID.video);
        room_GT(object_num(ii)) = objectID(ind,REF_ID.room);
        locInRoom_GT(object_num(ii)) = objectID(ind,REF_ID.locationInRoom);
        
        LocInRoomResponse_all(object_num(ii), subj_i) = LocInRoomResponse;
        
        % add post-hoc testing in the same form
        postHoc_response_hc(object_num(ii), subj_i) = objectID(ind,REF_ID.HC);
        postHoc_response_rc(object_num(ii), subj_i) = objectID(ind,REF_ID.RC);
        
    end

    allSubjectRating.Loc(subj) = nansum(objectID(LearnedSet,REF_ID.LocInRoomRespBallPark))/length(objid_str);
    allSubjectRating.ExactLocItemNum{subj} = []; 
    IDS = objectID(:,REF_ID.LocInRoomRespBallPark) == 1; % get rid of NaNs
    allSubjectRating.BallParkLocItemNum{subj} = objectID(IDS,REF_ID.object_num);
    
    % if sum(~ismember(allSubjectRating.ExactLocItemNum{subj} , allSubjectRating.BallParkLocItemNum{subj}))
    %    error('bad extraction of item numbers'); end
    allSubjectRating.LocInRoomRespBallPark(subj) = nansum(objectID(LearnedSet,REF_ID.LocInRoomRespBallPark))/length(objid_str);
    allSubjectRating.Rhit(subj) = nansum(objectID(LearnedSet,REF_ID.Rhit))/length(objid_str);
    allSubjectRating.HC(subj) = nansum(objectID(LearnedSet,REF_ID.HC))/length(objid_str);
    allSubjectRating.RC(subj) = nansum(objectID(LearnedSet,REF_ID.RC))/length(objid_str);
    
    allSubjectRating.RC_encodingCorrect(subj) = nansum( (objectID(LearnedSet,REF_ID.RC)== 1) & (objectID(LearnedSet,REF_ID.LocInRoomRespBallPark)== 1))/sum((objectID(LearnedSet,REF_ID.LocInRoomRespBallPark)== 1));
    
    % find items near front door
    IDS = find(objectID(LearnedSet,REF_ID.Rhit) == 1);
    allSubjectRating.rhit_items{subj} = objectID(LearnedSet(IDS),REF_ID.object_num);
   
    IDS = find(ismember(objectID(LearnedSet,REF_ID.locationInRoom),[1,2,5,6]) & ismember(objectID(LearnedSet,REF_ID.room),1));
    allSubjectRating.itemsNearFrontDoor{subj} = objectID(LearnedSet(IDS),REF_ID.object_num);
    IDS = find( (  ismember(objectID(LearnedSet,REF_ID.locationInRoom),7) &...
                    ismember(objectID(LearnedSet,REF_ID.room),1) &...
                    ismember(objectID(LearnedSet,REF_ID.house),[1,2])) | ...
                    (ismember(objectID(LearnedSet,REF_ID.locationInRoom),3) &...
                    ismember(objectID(LearnedSet,REF_ID.room),2) &...
                    ismember(objectID(LearnedSet,REF_ID.house),1))| ...
                    (ismember(objectID(LearnedSet,REF_ID.locationInRoom),2) &...
                    ismember(objectID(LearnedSet,REF_ID.room),2) &...
                    ismember(objectID(LearnedSet,REF_ID.house),2)) );            
    allSubjectRating.itemsNearInternalDoor{subj}  = objectID(LearnedSet(IDS),REF_ID.object_num);

    v1 = objectID(LearnedSet,REF_ID.RC);
    v1(isnan(v1)) = 0;
    v2 = objectID(LearnedSet,REF_ID.LocInRoomRespBallPark);
    v2(isnan(v2)) = 0;
    [RHO, P] = corr(v1,v2);
    allSubjectRating.RC_roomEncoding_RHO(subj) = RHO;
    
    subjResp_File = fullfile(folderName,sprintf('stimuliRespList_s%03d.mat',subj));
    save(subjResp_File, 'REF_ID','objectID')
    
    if PLOT
        figure(f1)
        cnt = cnt + 1;
        subplot(5,5,cnt);
        LocInRoomRespBallPark = objectID(LearnedSet,REF_ID.LocInRoomRespBallPark);
        LocInRoomRespBallPark(isnan(LocInRoomRespBallPark)) = 0;
        LocInRoomRespBallPark = logical(LocInRoomRespBallPark);
        remLoc = objectID(logical(LocInRoomRespBallPark),REF_ID.locationInRoom);
        N1 = nansum( LocInRoomRespBallPark & (objectID(LearnedSet,REF_ID.house)== 1 ));
        N2 = nansum( LocInRoomRespBallPark & (objectID(LearnedSet,REF_ID.house)==-1) );
        [N,C] = hist(remLoc,[1:8]);
        bar(C,N)
        set(gca,'xtick',[1,4,8])
        set(gca,'fontsize',7)
        title(sprintf('subj %d (rhit = %1.1f, room-correct %1.2f)',subj,allSubjectRating.Rhit(subj), length(remLoc)/200))
        
        figure(f2)
        subplot(5,5,cnt);
        LocInRoomRespBallPark = objectID(LearnedSet,REF_ID.LocInRoomResp) == objectID(LearnedSet,REF_ID.locationInRoom);
        LocInRoomRespBallPark(isnan(LocInRoomRespBallPark)) = 0;
        LocInRoomRespBallPark = logical(LocInRoomRespBallPark);
        remLoc = objectID(logical(LocInRoomRespBallPark),REF_ID.locationInRoom);
        N1 = nansum( LocInRoomRespBallPark & (objectID(LearnedSet,REF_ID.house)== 1 ));
        N2 = nansum( LocInRoomRespBallPark & (objectID(LearnedSet,REF_ID.house)==-1) );
        [N,C] = hist(remLoc,[1:8]);
        bar(C,N)
        hold all
        set(gca,'xtick',[1,4,8])
        set(gca,'fontsize',7)
        plot(get(gca,'xlim'),ones(1,2)*(1/8)*(200/8),'b--');
        title(sprintf('subj %d (rhit = %1.1f, sector-correct %1.2f)',subj,allSubjectRating.Rhit(subj), length(remLoc)/200))
        
    end
end
allSubjectRating.LocInRoomResponse_all = LocInRoomResponse_all;
allSubjectRating.postHoc_response_hc = postHoc_response_hc;
allSubjectRating.postHoc_response_rc = postHoc_response_rc;
    
itemNumber = 1:304;
items_info = table(itemNumber,houseType,videoNum, locInRoom_GT,room_GT);

% extract the items that were learned (non-lures)
learnedInd = find(sum(subj_view,2) > 0);
disp(length(learnedInd)) % should be 240 items


save(fullfile(onsetFiles_root,'subjectsEncodingRating.mat'),'allSubjectRating',...
        'items_info','learnedInd','subj_view', 'subj_correct', 'subj_incorrect', 'houseType', 'room_GT')

% remove NaNs

if PLOT
    newA4figure('behavStats')
    bins = [0:0.05:1];
    subplot(2,2,1)
    hold all
    histogram(allSubjectRating.Rhit(subjects),bins)
    title('rhit in scanner')
    plot(ones(1,2)*1/2, get(gca,'ylim'),'b--')
    subplot(2,2,2)
    hold all
    histogram(allSubjectRating.HC(subjects),bins)
    plot(ones(1,2)*1/2, get(gca,'ylim'),'b--')
    title('House Correct - post scan')
    subplot(2,2,3)
    hold all
    hold all
    histogram(allSubjectRating.RC(subjects),bins)
    plot(ones(1,2)*1/2, get(gca,'ylim'),'b--')
    title('Room Correct - post scan')
    
    subplot(2,2,4)
    hold all
    histogram(allSubjectRating.LocInRoomRespBallPark(subjects),bins)
    plot(ones(1,2)*1/8, get(gca,'ylim'),'b--')
    title('Room Correct - encoding phase')
    
    
    [rho(1),pval(1)] = corr(allSubjectRating.Rhit(subjects)',allSubjectRating.HC(subjects)','Type','Spearman');
    [rho(2),pval(2)] = corr(allSubjectRating.Rhit(subjects)',allSubjectRating.Loc(subjects)','Type','Spearman');
    [rho(3),pval(3)] = corr(allSubjectRating.HC(subjects)',allSubjectRating.Loc(subjects)','Type','Spearman');
    [rho(4),pval(4)] = corr(allSubjectRating.Rhit(subjects)',allSubjectRating.RC(subjects)','Type','Spearman');
    [rho(5),pval(5)] = corr(allSubjectRating.LocInRoomRespBallPark(subjects)',allSubjectRating.RC(subjects)','Type','Spearman');
    
    title('Location pretty accurate')
    
    figure; hold all
    pp = plot(allSubjectRating.LocInRoomRespBallPark(subjects)',allSubjectRating.Rhit(subjects)','.');
    pp.Color = 'k';
    pp.MarkerSize = 15;
    axis([0.8 1 0 1])
    plot(ones(1,2)*1/8,get(gca,'ylim'),'b-');
    plot(get(gca,'xlim'),ones(1,2)*.5,'b-');
    xlabel('Prob for accurate sector during encoding')
    ylabel('Prob for recognition in the final test')
    legend('subject')
    
    
    figure; hold all
    pp = plot(allSubjectRating.LocInRoomRespBallPark(subjects)',allSubjectRating.RC(subjects)','.');
    pp.Color = 'k';
    pp.MarkerSize = 15;
    axis([0.8 1 0 1])
    plot(ones(1,2)*0.5,get(gca,'ylim'),'b-');
    plot(get(gca,'xlim'),ones(1,2)*.5,'b-');
    xlabel('Prob for accurate sector during encoding')
    ylabel('Prob for correct room in the final test')
    legend('subject')
    title(sprintf('RHO = %2.2f, p = %2.2e',rho(5),pval(5)))
    
    figure; hold all
    pp = plot(allSubjectRating.LocInRoomRespBallPark(subjects)',allSubjectRating.RC_encodingCorrect(subjects)','.');
    hold all
    pp.Color = 'k';
    pp.MarkerSize = 15;
    axis([0.7 1 0 1])
    plot(ones(1,2)*.5,get(gca,'ylim'),'b-');
    plot(get(gca,'xlim'),ones(1,2)*.5,'b-');
    ylabel('(Encoding & post successful) /(Encoding successful)')
    xlabel('Prob for accurate room during encoding')
    legend('subject')
    
    
end

PrintActiveFigs(outputFigures)
end