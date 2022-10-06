% Adjustment of code from (intended for use with DML_batch_qa) to compare
% SNR for additional sub-fields and MTL areas used in this project.

function [b,snrmean, snrmean_mid, dvars, spatialsnr_mean, GS_suspects]=compute_snr_spatcon_ROIs(b,figFormat,pdfReport,GSthresh)
% This function generates SNR maps and generates and saves plots of
% intensity over time and mean SNR by run
% extended from scripts from here: http://dbic.dartmouth.edu/wiki/index.php/Noise_Detection
%
% Author: Maureen Ritchey, Shao Fang Wang, Liang-Tien (Frank) Hsieh,
% July 2014
% Intended for use with DML_batch_qa

b.QA_dir_root = fullfile(b.dataDir, sprintf('QA_main_runs'));
subj_roi_dir = fullfile('E:\spatcon_hdz\data_p\ROIs_from_OSF',sprintf('%s_rois',b.curSubj),sprintf('%s',b.curSubj),'ROIs');
roi_title = {'ca1','ca23dg','erc','sub','alerc','pmerc','phcant'};
roi_list = {fullfile(subj_roi_dir,'ashs_left','brCA1_body.nii'),...
    fullfile(subj_roi_dir,'ashs_left','brCA2_3_DG_body.nii'),...
    fullfile(subj_roi_dir,'ashs_left','brERC.nii'),...
    fullfile(subj_roi_dir,'ashs_left','brsubiculum.nii'),...
    fullfile(subj_roi_dir,'manual','bralERC_L.nii'),...
    fullfile(subj_roi_dir,'manual','brpmERC_L.nii'),...
    fullfile(subj_roi_dir,'manual','brPHC_Ant_L.nii')};
    
%Load brain mask
whole_brain_mask_filename = b.mask;
M = spm_vol(b.mask);
whole_brain_mask = spm_read_vols(M);

for ii = 1:length(roi_list)
    
    b.QA_dir = fullfile(b.QA_dir_root,roi_title{ii});
    mkdir(b.QA_dir);
    b.mask = roi_list{ii};
    
    if isempty(dir(b.mask))
        disp(b.mask)
        disp('missing mask')
        continue
    end
    %Load brain mask
    M = spm_vol(b.mask);
    mask = spm_read_vols(M);
    
    
    %% INITIALIZE VARIABLES ETC
    snrmean = []; snrmean_mid = [];
    close all;
    f = figure('Visible','off');
    st = figure('Visible','off');
    g1 = figure('Visible','off');
    g2 = figure('Visible','off');
    g3 = figure('Visible','off');
    s = figure('Visible','off');
    dv = figure('Visible','off');
    gs = figure('Visible','off');
    sp = figure('Visible','off');
    
    %% LOOP OVER RUNS
    % Start loop through runs
    
    
    
    %Initialize global signal suspects
    GS_suspects = {};
    
    %Loop over runs
    for j=1:length(b.runs)
        fprintf('\n%s (%02d / %02d runs)\n',b.runs{j},j,length(b.runs))
        
        %% DATA HANDLING
        
        %Load preprocessed BOLDS
        disp('Preparing images')
        P = b.rundir(j).rfiles;
        files = spm_vol(P);
        
        %Initialize variables
        outputname = b.runs{j};% current run name
        data = spm_read_vols(files(1));
        avg= zeros(size(data));% create an empty variable sized by the image dimensions
        sd_tmp=zeros(size(data));% again, creating empty variable sized by the image dimensions
        ntimepoints = size(files,1);% figure out how many files you're dealing w
        slicematrix = [];
        
        %Store data and calculate metrics
        disp('Loading images')
        for i = 1:ntimepoints % loop across .nii files
            % concatenate all images into a 4 D matrix. The 4th dimension is time (n).
            data(:,:,:,i) = spm_read_vols(files(i));
        end
        avg=sum(data,4)/ntimepoints; % average across time points
        
        %% TEMPORAL SNR
        disp('Computing temporal SNR')
        for i=1:ntimepoints
            sd_tmp(:,:,:,i)=(avg-double(data(:,:,:,i))).^2;
            % average slices over within-brain voxels
            slicematrix=[ slicematrix squeeze(mean(mean(mask.*data(:,:,:,i))))];
        end
        sd_tmp=sum(sd_tmp,4);
        
        sd=sqrt(sd_tmp/(ntimepoints-1));% standard SD calculation
        snr=avg./sd;% standard SNR calculation
        meantimecourse = mean(slicematrix);
        
        %Clean up SNR variable
        snr = mask.*snr; % mask out non-brain data
        snr(snr>2000) = NaN; % mask out abnormal values
        snr(snr==0) = NaN;
        snrmean = [snrmean nanmean(nanmean(nanmean(snr(:,:,:))))];
        snrmean_mid = [snrmean_mid nanmean(nanmean(snr(:,:,(size(data,3)/2))))];
        snr(isnan(snr))=0;
        
        %% CALCULATE FRAME-WISE SIGNAL STATS
        tmpdata = reshape(data,numel(data(:,:,:,1)),size(data,4))';
        tmpmask = reshape(mask,numel(data(:,:,:,1)),1)';
        tmpdata(:,tmpmask==0 | sum(tmpdata)==0) = []; %remove voxels outside the brain
        
        % Calculate DVARS (as in Power et al. 2012)
        % Equation: DVARS(i) = sqrt(brainaverage((intensity(i) - intensity(i-1))^2)))
        disp('Computing DVARS')
        % Convert to mode1000 normalized data
        modeval = mode(reshape(tmpdata,numel(tmpdata),1)); % get mode across all voxels and timepoints
        tmpdata_mn = (1000/modeval).*tmpdata; % normalize: mode=1000, 10-point change=1%
        deltaIntensity = [zeros(1,size(tmpdata_mn,2)); diff(tmpdata_mn)];
        % note that the next 2 steps can be equivalently replaced with rms(deltaIntensity')'
        avgDeltaIntensity = nanmean(deltaIntensity.*deltaIntensity,2);
        dvars = sqrt(avgDeltaIntensity);
        
        % Calculate percent change from global mean (as in ART repair toolbox)
        disp('Computing percent change from global mean')
        meanval = mean(reshape(tmpdata,numel(tmpdata),1));
        tmpdata_g = (100/meanval).*tmpdata;
        globalsignal = nanmean(tmpdata_g,2) - 100;
        clear tmpdata tmpmask tmpdata_mn tmpdata_g deltaIntensity deltaIntensity2 avgDeltaIntensity
        
        % Identify bad timepoints based on global mean intensity
        GS_suspects{j} = find(abs(globalsignal)>GSthresh);
        fprintf('%0.0f bad timepoints (%0.0f percent) identified based on global signal\n\n', ...
            length(GS_suspects{j}),100.*length(GS_suspects{j})./size(globalsignal,1));
        
        %output files to the QA dir
        outputdir = b.QA_dir;
        avg_output=files(1);
        sd_output=files(1);
        snr_output=files(1);
        avg_output.fname    = [outputdir,outputname,'_average.nii'];
        sd_output.fname    = [outputdir,outputname,'_SD.nii'];
        snr_output.fname    = [outputdir,outputname,'_SNR.nii'];
        %     spm_write_vol(avg_output,avg); % typically don't need these but can uncomment if you want them
        %     spm_write_vol(sd_output,sd);
        spm_write_vol(snr_output,snr);
        
        %% SPATIAL SNR
        % Compute spatial SNR
        fprintf('Computing spatial SNR\n')
        clear signal_mean noise_SD;
        tmpdata = reshape(data,numel(data(:,:,:,1)),size(data,4))';
        tmpmask = reshape(mask,numel(data(:,:,:,1)),1)';
        tmpmask_whole_brain = reshape(whole_brain_mask,numel(data(:,:,:,1)),1)';
        
        signaldata = tmpdata; noisedata = tmpdata;
        signaldata(:,tmpmask==0 | sum(tmpdata)==0) = []; %remove voxels outside the requested ROI
        noisedata(:,tmpmask_whole_brain==1 | sum(tmpdata)==0) = []; %keep voxels outside the brain
        signal_mean = mean(signaldata'); %will return a row vector w/ the number of elements == number of images
        noise_SD    = std(double(noisedata')); %change to a double variable type (for increased precision) and take SD
        spatial_SNR = 0.655*(signal_mean./noise_SD);%constant (from the web) for fmri. then divide every element of the vector by SD_noise
        % The factor of 0.655 arises because magnitude images have a Rician noise distribu&on which
        % has no negative values (unlike Gaussian noise for complex signals)
        spatial_SNR_runs{j}=spatial_SNR;
        clear tmpdata tmpmask signaldata noisedata
        
        %% MAKE PLOTS
        %Save jpgs of slices
        avg_slices = save_slices(avg,[outputdir,outputname,'_average.jpg'],0);
        sd_slices = save_slices(sd,[outputdir,outputname,'_SD.jpg'],0);
        snr_slices = save_slices(snr,[outputdir,outputname,'_tSNR.jpg'],0);
        
        %plot average signal for middle slice over time
        disp('Plotting mean time course, slice timepoint matrix, average slices, sd slices, and snr slices')
        figure(f)
        subplot(length(b.runs),1,j);
        plot(meantimecourse)
        line([1 size(meantimecourse,2)],[mean(meantimecourse)+3.*std(meantimecourse) mean(meantimecourse)+3.*std(meantimecourse)],'Color','k');
        line([1 size(meantimecourse,2)],[mean(meantimecourse)-3.*std(meantimecourse) mean(meantimecourse)-3.*std(meantimecourse)],'Color','k');
        axis([0 size(meantimecourse,2) -2 2]); axis 'auto y';
        xlabel('timepoints'); ylabel('au');
        title(['Mean Intensity: ',b.curSubj,' ',b.runs{j}],'Interpreter','none');
        
        %plot slice x timepoint matrix
        figure(st)
        subplot(length(b.runs),1,j);
        imagesc(slicematrix)
        xlabel('timepoints'); ylabel('slices');
        title(['Slice x Timepoint: ',b.curSubj,' ',b.runs{j}],'Interpreter','none');
        colorbar;
        
        %plot average slices
        figure(g1)
        subplot(length(b.runs),1,j);
        image(avg_slices)
        colormap(gray)
        axis off
        title(['Mean: ',b.curSubj,' ',b.runs{j}],'Interpreter','none');
        
        %plot sd slices
        figure(g2)
        subplot(length(b.runs),1,j);
        image(sd_slices)
        colormap(jet)
        axis off
        title(['SD: ',b.curSubj,' ',b.runs{j}],'Interpreter','none');
        
        %plot snr slices
        figure(g3)
        subplot(length(b.runs),1,j);
        image(snr_slices)
        colormap(gray)
        axis off
        title(['tSNR: ',b.curSubj,' ',b.runs{j}],'Interpreter','none');
        
        %plot dvars over time
        figure(dv)
        subplot(length(b.runs),1,j);
        plot(dvars)
        axis([1 size(dvars,1) -2 2]); axis 'auto y';
        xlabel('timepoints'); ylabel('DVARS');
        title(['DVARS: ',b.curSubj,' ',b.runs{j}],'Interpreter','none');
        
        %plot global signal over time
        figure(gs)
        subplot(length(b.runs),1,j);
        plot(globalsignal)
        line([1 size(globalsignal,1)],[GSthresh GSthresh],'Color','r');
        line([1 size(globalsignal,1)],[-GSthresh -GSthresh],'Color','r');
        axis([1 size(globalsignal,1) -2 2]); axis 'auto y';
        xlabel('timepoints'); ylabel('GS %');
        title(['Global signal (% from mean): ',b.curSubj,' ',b.runs{j}],'Interpreter','none');
        
        clear data avg snr sd files
    end
    
    %% SAVE PLOTS
    
    %finish timecourse plots
    disp('Saving plots')
    save_figs(f,'timecourse_mean_plot',b,figFormat,pdfReport);
    save_figs(st,'timecourse_slice_plot',b,figFormat,pdfReport);
    save_figs(g1,'slices_avg_plot',b,figFormat,pdfReport);
    save_figs(g2,'slices_sd_plot',b,figFormat,pdfReport);
    save_figs(g3,'slices_tsnr_plot',b,figFormat,pdfReport);
    save_figs(dv,'dvars_plot',b,figFormat,pdfReport);
    save_figs(gs,'global_signal_plot',b,figFormat,pdfReport);
    
    %plot mean SNR by run
    figure(s)
    snrmean = [snrmean' snrmean_mid'];
    bar(snrmean);
    title([b.curSubj '-Temporal SNR'],'FontSize',14);
    xlabel('run'); ylabel('mean SNR');
    legend('whole-brain','midslice','Location','EastOutside');
    save_figs(s,'tsnr_plot',b,figFormat,pdfReport,6,6); % square plot
    
    % % Plot spatial SNR
    figure(sp)
    maxrunlen = size(spatial_SNR,2);
    colors = {'r','g','b','k','y','m','c'};
    for j=1:size(spatial_SNR_runs,2) % coded this way in case runs have different lengths
        if size(spatial_SNR_runs{j},2) > maxrunlen
            maxrunlen = size(spatial_SNR_runs{j},2);
        end
        plot(spatial_SNR_runs{j},colors{(mod(j-1,length(colors))+1)}); % cycle through colors
        hold on
    end
    hold off
    legend(b.runs,'Location','EastOutside')
    title([b.curSubj '-Spatial SNR'],'FontSize',14);
    ylabel('Spatial SNR','FontSize',14)
    xlabel('Timepoint','FontSize',14)
    axis([0 maxrunlen -2 2]); axis 'auto y';
    save_figs(sp,'spatialSNR',b,figFormat,pdfReport,6,6); % square plot
    
    %Compute mean spatial_SNR for summary info
    for j=1:size(spatial_SNR_runs,2)
        spatialsnr_mean(j) = mean(spatial_SNR_runs{j}); % mean(spatial_SNR_runs,2);
    end
    
    filename = fullfile(b.QA_dir,sprintf('QA_result_%s_left_%s',b.curSubj,roi_title{ii}));
    save(filename,'spatialsnr_mean','snrmean','dvars')
    b.QA_dir = b.QA_dir_root;
     
end %snr_temporal

b.mask = whole_brain_mask_filename;


end

%% SUBFUNCTIONS
function [allslices]=save_slices(data,outputFile,writeFlag)
%this function takes a 3-D image and horizontally concatenates the select
%slices so they can be viewed in a row; returns the data for use in a plot,
%and optionally saves the multi-slice image
allslices = [];
for j=1:3:size(data,3) % shows every third slice
    slice = data(:,:,j)';
    allmin = nanmin(nanmin(nanmin(data(:,:,:))));
    allmax = nanmax(nanmax(nanmax(data(:,:,:))));
    slice=round(63*(slice-allmin)/allmax)+1;
    allslices = [allslices slice];
end
if writeFlag,imwrite(allslices,gray,outputFile);end

end %save_slices