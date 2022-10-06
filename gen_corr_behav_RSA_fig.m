% Create figure that demonstrates a correlation between behavior and fMRI
% corr analysis
% Depends on gen_corr_behav_RSA_fig.m 

function gen_corr_behav_RSA_fig()

% prep data
mixed_models_add_behav_for_fig()

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

%%
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

%%
newA4figure('RSA_conservative_diff_corr_behav_rhit')
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
C = RSA_score_SR_SUBJ_ROI(ind,hem_i,roi_i);

hold all
for ii = 1:length(rhit_perc_all)
    plot(C(ii),100*B(ii),'o','color',ccmap(5+ii,:),'linewidth',2)
    if pval_subj_v2(ii,hem_i,roi_i)<0.05
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
   sprintf('%s, R = %2.2f, p = %2.2e',strrep(roi_str{roi_i},'_',' '), R_rhit_SR(hem_i,roi_i),corr_p_rhit_SR(hem_i,roi_i)),...
   'HorizontalAlignment','center');
xlabel('RSA spatial sensitivity score')
ylabel('Recollection (%)')

coefs = polyfit(C, 100*B', 1);
% Plot line
ah = refline(coefs); 
ah.Color = 'k';
set(gca,'ylim',100*[0 1])

end

%%
PrintActiveFigs('E:\Dropbox\RanganathLab\spatcon_2\figures\LINKS');


end