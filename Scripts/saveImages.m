function saveImages(saveDir,saveName,imgTrace)
if(~exist(saveDir, 'dir'))
    mkdir(saveDir)
end
if(strcmp(get(gcf,'type'),'figure'))
    set(gcf,'Renderer','painters');
end
 axis image               % resolution based on image
 axis off                 % avoid printing axis
 set(gcf, 'color', 'none');
 set(gca, 'color', 'none');
 set(gca,'LooseInset',get(gca,'TightInset')); % removing extra white space in figure
 sName = [saveDir,saveName,'.eps'];
exportgraphics(gcf,sName,'ContentType','vector');
if(~isempty(imgTrace))
    boundaries = bwboundaries(imgTrace);
    fHandle2 = copyobj(gca(fHandle),figure()); hold on;
    cellfun(@(b) plot(b(:,2), b(:,1), 'g', 'LineWidth',lineWidth), boundaries);
    saveFigures(fHandle2,[saveDir,'Imaging\'],saveName,[]);
end
end