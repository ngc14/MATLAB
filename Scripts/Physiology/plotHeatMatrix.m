
function plotHeatMatrix(hm,hmNames,entryNames)
colormap('jet');
for c = 1:length(hm)
    nexttile; hold on; axis image; axis ij;
    title(hmNames(c));
    imagesc(hm{c});
    xticks(1:length(entryNames));xticklabels(entryNames);yticklabels([]);clim([0 1]);
end
end
% subplot(2,length(conditions)+2,(length(conditions)+2)+1);hold on; axis image; axis ij;
% title('Go correlations');
% imagesc(goMatrix);
% xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
% subplot(2,length(conditions)+2,(length(conditions)+2)+4); hold on; axis ij; axis tight;
% title('Go marginals');
% imagesc(mean(goMatrix.*(~logical(diag(true(sum(taskUnits),1)))./~logical(diag(true(sum(taskUnits),1)))),2,'omitnan'));
% xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);clim([0 .5]);
% subplot(2,length(conditions)+2,(length(conditions)+2)+2);hold on; axis image; axis ij;
% title('Reach correlations');
% imagesc(reachMatrix);
% xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
% subplot(2,length(conditions)+2,(length(conditions)+2)+5); hold on; axis ij; axis tight;
% title('Reach marginals');
% imagesc(mean(reachMatrix.*(~logical(diag(true(sum(taskUnits),1)))./~logical(diag(true(sum(taskUnits),1)))),2,'omitnan'));
% xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);clim([0 .5]);
% subplot(2,length(conditions)+2,(length(conditions)+2)+3); hold on; axis ij; axis image;
% title('Grasp correlations');
% imagesc(graspMatrix);
% xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
% subplot(2,length(conditions)+2,(length(conditions)+2)+6); hold on; axis ij; axis tight;
% title('Grasp marginals');
% imagesc(mean(graspMatrix.*(~logical(diag(true(sum(taskUnits),1)))./~logical(diag(true(sum(taskUnits),1)))),2,'omitnan'));
% xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);colorbar;clim([0 .5]);
%
%
%
% 
% subplot(1,length(conditions)+3,c);hold on; axis image; axis ij;
% title(conditions(c));
% imagesc(concatMatrix{c});
% numPairs{c} = sum(~isnan(condAllMatrix),3);
% xticks(1:2:32);xticklabels(arrayfun(@num2str,1:2:32,'UniformOutput',false));yticklabels([]);clim([0 1]);
% avgMatrix = mean(cat(3,concatMatrix{1:length(conditions)-1}),3,'omitnan');
% as=subplot(1,length(conditions)+3,length(conditions)+1); hold on; axis image; axis ij;
% title("Average");
% imagesc(avgMatrix);
% xticks(1:2:32);xticklabels(arrayfun(@num2str,1:2:32,'UniformOutput',false));yticklabels([]);clim([0 1]);
% subplot(1,length(conditions)+3,length(conditions)+3); hold on; axis ij; axis tight;
% imagesc(mean(avgMatrix.*(~logical(diag(1:32))./~logical(diag(1:32))),2,'omitnan'));
% title('Average marginals');
% colorbar;clim([0 1]);xlim([0.5 1]);xticklabels([]);yticks(1:32);yticklabels(arrayfun(@num2str,1:32,'UniformOutput',false));
