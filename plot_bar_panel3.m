function P_mat = plot_bar_panel3(Z,xlabel_str,barcolor)

if nargin < 3 % default is blue bars
    barcolor(1,:) = [0,0,1];
    barcolor(2,:) = [0,0,1];
end

if length(xlabel_str) > 3; error('too many tick labels'); end
if sum(size(Z) ~= [1,3]); error('Z size off'); end

for ii = 1:3
    meanVal(ii) = nanmean(Z{ii});
    semVal(ii) = nanstd(Z{ii})/sqrt(length(Z{ii}));
    bb = bar(ii,meanVal(ii),'FaceColor',barcolor(ii,:));
    hold all
end
errorbar(1:3,[meanVal(1),meanVal(2),meanVal(3)],...
            [semVal(1),semVal(2),semVal(3)],'k.')
        
[P_mat(1), h] = ranksum(Z{1},Z{2});
[P_mat(2), h] = ranksum(Z{2},Z{3});
[P_mat(3), h] = ranksum(Z{1},Z{3});
[P_mat(4), h] = ranksum(Z{1},[Z{3}',Z{2}']);


box off

% if P_mat(1) < 0.05
%     plot(1.5,max(medianVal(:))*1.07,'k*'); end

set(gca,'XTick',1:3,'XTickLabel',xlabel_str,'xticklabelrotation',35)

end