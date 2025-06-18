function saveFigures(fHandle,saveDir,saveName,imgTrace)
if(~exist(saveDir, 'dir'))
    mkdir(saveDir)
end
if(strcmp(get(fHandle,'type'),'figure'))
    set(fHandle,'Renderer','painters');
end
exportgraphics(fHandle,strcat(saveDir,saveName,'.eps'),'ContentType','vector');
pause(1);
exportgraphics(fHandle,strcat(saveDir,saveName,'.png'));
pause(1);
savefig(fHandle,strcat(saveDir,saveName));
%saveas(fHandle,[saveDir,saveName,'.fig']);
if(~isempty(imgTrace))
    boundaries = bwboundaries(imgTrace);
    fHandle2 = copyobj(gca(fHandle),figure()); hold on;
    cellfun(@(b) plot(b(:,2), b(:,1), 'g', 'LineWidth',lineWidth), boundaries);
    saveFigures(fHandle2,[saveDir,'Imaging\'],saveName,[],'ContentType','vector');
end
end