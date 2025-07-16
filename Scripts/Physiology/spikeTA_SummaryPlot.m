monkey ='Gilligan';
sessionDate ='06_13_2019';
smoothKernel = 5;
debug = false;
plotISA = true;
testWindow = [0];
preWindow = [-20 -5];
if(strcmp(monkey, 'Gilligan'))
    sessionDateF = string(datetime(sessionDate,'InputFormat','MM_dd_yyyy','Format','MM_dd_yyyy'));
else
    sessionDateF = string(datetime(sessionDate,'InputFormat','yyyy_MM_dd','Format','yyyy_MM_dd'));
end
saveDir = 'S:\Lab\ngc14\Working\EMG_UNITS\Summarys\';
load(strcat(saveDir,sessionDateF,"\unitXmuscle.mat"));
Fs = 2000;
muscles = [{'Deltoid'}, {'Biceps'},{'Triceps'},{'Wrist extensor'},{'Wrist flexor'},{'Digit extensor'},{'Digit flexor'}];
cmap = [];
cmap(1,:) = [.7 0 .8];
cmap(2,:) = [.9 .97 .9];
cmap(3,:) = [1 .7 0];
[X,Y] = meshgrid(1:3,1:255);
cmap = interp2(X([1,128,255],:),Y([1,128,255],:),cmap,X,Y);

[~,~,testInd] = intersect(testWindow./(1000/Fs),indexWindow(1):indexWindow(end));
[~,~,preInd] = intersect(preWindow./(1000/Fs),indexWindow(1):indexWindow(end));
subFrags = cellfun(@(s) s{1}-s{2}, unitXmuscle, 'UniformOutput',false);
difs = cellfun(@(s) [{mean(s,1,'omitnan')},{mean(s(:,preInd(1):preInd(end)),2,'omitnan')}], ...
    subFrags,'UniformOutput',false);
aboveSTD = (cellfun(@(s) find(mean(s{2},'omitnan') + 2*std(s{2},'omitmissing') < s{1}),difs,'UniformOutput',false));
belowSTD = (cellfun(@(s) find(mean(s{2},'omitnan') - 2*std(s{2},'omitmissing') > s{1}),difs,'UniformOutput',false));
[rSig,cSig]= find(cellfun(@any,aboveSTD)|cellfun(@any,belowSTD));
boxDists = cellfun(@(m) [m{1},m{2}],unitXmuscle, 'UniformOutput',false);
summary = NaN(2,size(difs,2));
[~,indMat] = cellfun(@(s) max(abs(s{1}(:,testInd(1):end))),difs,'UniformOutput',false);
testMat = cell2mat(cellfun(@(s,b) reshape(s{1}(:,testInd(1)+b-1),1,1,[]),difs,indMat, 'UniformOutput',false));
baseMat = cell2mat(cellfun(@(s) reshape(mean(s{2},'omitnan'),1,1,[]),difs,'UniformOutput',false));
sigVals = (testMat-baseMat).*double(cellfun(@any,aboveSTD)|cellfun(@any,belowSTD));
summary(1,:) = mean(sigVals.*double((sigVals>0))./double(sigVals>0),1,'omitnan');
summary(2,:) = mean(sigVals.*double((sigVals<0))./double(sigVals<0),1,'omitnan');
cl = [ median(cell2mat(cellfun(@(s) reshape(mean(s{2},'omitnan')-(2*std(s{2},'omitmissing')),1,1,[]),difs,'UniformOutput',false)),'all'),...
     median(cell2mat(cellfun(@(s) reshape(mean(s{2},'omitnan')+(2*std(s{2},'omitmissing')),1,1,[]),difs,'UniformOutput',false)),'all')];
f4 = figure(4);
imagesc([testMat-baseMat;-Inf(1,size(difs,2));summary]);
hold on;
axis IJ;
colormap([0 0 0 ;cmap]);
clim(cl);
xticks(1:length(muscles));
xticklabels(muscles);
yticks(1:size(difs,1)+3);
yticklabels([vertcat(reshape(cellstr(arrayfun(@(s,u) string(strcat(num2str(s),",",num2str(u))), spkChannel,unitNum,'UniformOutput',false)),[],1),...
    "","Facilitation","Suppression")])
title('Max Post-Spike');
colorbar;
imH = imhandles(f4);
for s = 1:length(rSig)
    plot(gca(f4),[cSig(s)-.5,cSig(s)-.5,cSig(s)+.5,cSig(s)+.5,cSig(s)-.5],...
        [rSig(s)-.5,rSig(s)+.5,rSig(s)+.5,rSig(s)-.5,rSig(s)-.5],'k-','LineWidth',2);
end

saveFigures(f4,strcat(saveDir,sessionDateF,'\'),"Summary_MAX",[]);
%% debug distributions
if(debug)
    f5=gca(figure(5));
    while 1
        ch = drawpoint(imH.Parent);
        pos = round(ch.Position);
        col = pos(1);
        row = pos(2);
        boxplot((f5),boxDists{row,col}(:,[1 2]),'Notch','on');
        xticklabels((f5),{'Test','Baseline'});
        ch.Visible = 'off';
    end
end
