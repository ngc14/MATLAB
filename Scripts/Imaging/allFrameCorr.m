monkey = 'Gilligan';
conds = {'[ExtraSmallSphere]','[LargeSphere]', '[Photocell]','[Rest]'};
if(strcmp(monkey,'Gilligan'))
    mm = MotorMapping(42);
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
    frames = {49:53}
elseif(strcmp(monkey,'Skipper'))
    mm = MotorMapping(56);
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
    frames = {48:52}
end
HPk = 250;
LPk = 5;

condComb = nchoosek(1:length(conds),2);
condComNames = nchoosek(string(conds),2);
condComb(end+1:end+length(conds),:) = repmat(1:length(conds),[2 1])';
condComNames(end+1:end+length(conds),:) = repmat(string(conds),[2 1])';

[datesTable, masksR, ~] = getMonkeyInfo('S:\Lab\',monkey,["M1", "PMd"],false);
refMask = masksR{1};
siteMask = getVoronoiMask(datesTable,mm,refMask,["Arm","Hand"]);
monkeySiteMask = (siteMask & refMask) ./ siteMask;
%%
[evenCond, oddCond] = deal(cell(1,length(conds)));
tMat = cell(1,length(dates));
spacing = cell(1,length(dates));
for c = 1:length(conds)
    currCond = conds{c};
    eDates = {};
    oDates = {};
    for d = 1:length(dates)
         if (c==1)
            tform = matfile(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},...
                '\Imaging\run0',num2str(runs{d}),'\tform_nonrigid.mat']);
            Otrans{d} = tform.O_trans;
            spacing{d} = tform.Spacing;
            tform = [];
        end
    end
    for f = 1:length(frames)
        eFrame = {};
        oFrame ={};
        fT1 = datetime('now');
        parfor (d = 1:length(dates), floor(length(dates)/2))
            fPath = ['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},...
                '\Imaging\run0',num2str(runs{d}),'\Results\Blockview\',currCond,...
                '\FF_1\LP',num2str(LPk),'_HP',num2str(HPk),'\Align_Offset\noClip\'];
            if(exist([fPath, 'NaN_filtered_HP',num2str(HPk),'_LP',num2str(LPk),'.h5'],'file'))
                currTrials = squeeze(mean(h5read([fPath,'NaN_filtered_HP',...
                    num2str(HPk),'_LP',num2str(LPk),'.h5'],'/allTrials',...
                    [1 1 frames{f}(1) 1],[size(refMask,1),size(refMask,2),length(frames{f}),Inf]),3,'omitnan'));
                eFrame{d} = bspline_transform(Otrans{d},mean(currTrials(:,:,2:2:size(currTrials,3)),3,'omitnan'), spacing{d},3);
                oFrame{d} = bspline_transform(Otrans{d},mean(currTrials(:,:,1:2:size(currTrials,3)),3,'omitnan'), spacing{d},3);
            end
        end
        eDates{f} = eFrame;
        oDates{f} = oFrame;
    end
    evenCond{c} = eDates;
    oddCond{c} = oDates;
end
% for d = 1:size(runs,2)
%     sessionIms = cellfun(@(i,r) i(:,:,r{d}), ims,dateLabels, 'UniformOutput',false);
%     normR{d}{1} = nanmean(cat(3,sessionIms{:}),3);
%     normR{d}{2} = nanstd(cat(3,sessionIms{:}),0,3);
%     if(any(normR{d}{2}==0,'all'))
%         [r,c] = find(normR{d}{2}==0);
%         normR{d}{2}(r,c) = nanmin(abs(normR{d}{2}(normR{d}{2}~=0)));
%     end
%     normR{d}{1} = zeros(size(refMask));
%     normR{d}{2} = ones(size(refMask));
% end
%%
corrCondFrames = cell(1,length(condComb));
pairN = nchoosek(["e1","e2","o1","o2"],2);
condPairInds = ~(sum(contains(pairN, num2str(1)),2)==2 | sum(contains(pairN, num2str(2)),2)==2);
for c = 1:length(condComb)
    corrFrames = cell(1,length(frames));
    for f = 1:length(frames)
        if(condComb(c,2)==condComb(c,1))
            corrPairs = {oddCond{condComb(c,1)}{f},evenCond{condComb(c,2)}{f}};
        else
            c1E = evenCond{condComb(c,1)}(f);
            c2E = evenCond{condComb(c,2)}(f);
            c1O = oddCond{condComb(c,1)}(f);
            c2O = oddCond{condComb(c,2)}(f);
            corrPairs = nchoosek([c1E, c2E, c1O, c2O],2);
            corrPairs = corrPairs(condPairInds,:);
        end
        cMat = {};
        parfor m =1:size(corrPairs,1)
            c1 = corrPairs{m,1};
            c2 = corrPairs{m,2};
            c1N = cellfun(@isempty,c1);
            c2N = cellfun(@isempty,c2);
            c1(c1N) = {NaN(size(refMask))};
            c2(c2N) = {NaN(size(refMask))};
            c1 = cellfun(@(ee) ee.*monkeySiteMask, c1, 'UniformOutput',false);
            c2 = cellfun(@(ee) ee.*monkeySiteMask, c2, 'UniformOutput',false);
            cMat{m} = cell2mat(cellfun(@(cs1) cellfun(@(cs2) xcorr(cs1(~isnan(cs1)),...
                cs2(~isnan(cs2)),0,'normalized'),c2,'UniformOutput',true), c1,'UniformOutput',false)');
            cMat{m}(cMat{m}==0) = NaN;
%             cMat(c1N | c2N) = NaN;
        end
        corrFrames{f} = cat(3,cMat{:});

        %         if(sum(c2E)~=sum(c1E))
%             if(sum(c2E)>sum(c1E))
%                 c1E = c2E;
%             else
%                 c2E = c1E;
%             end
%         elseif(find(c1E)~=find(c2E))
%             c1E = c1E | c2E;
%             c2E = c1E | c2E;
%         end
    end
    corrCondFrames{c} = corrFrames;
end
%%
save(['S:\Lab\ngc14\Working\',monkey,'\XCorr\X.mat'],'corrCondFrames','-v7.3');
%%
corrCondFrames = cellfun(@(gp) cellfun(@(p) mean(p,3, 'omitnan'), gp, 'UniformOutput', false), corrCondFrames, 'UniformOutput',false);
corrCondFrames = cellfun(@(gp) cell2mat(cellfun(@(p) p(logical(eye(size(p,1)))), gp, 'UniformOutput',false))', corrCondFrames, 'UniformOutput',false);
groups = zeros(1,length(condComb));
groups(:,strcmp(condComNames(:,1), condComNames(:,2))) = 1;
groups(:,sum(strcmp(condComNames,"[Rest]"),2)==1) = 2;
groups = groups + 1;
colorsG = {'k','m','g','b','c','r','k','g','c',[1 .5 0],'r'};
%%
for n = 1:max(groups)+1
    ga = gca(figure());
    hold on;
    if(n<=max(groups))
        plotG = corrCondFrames(groups==min(n,max(groups)));
    else
        plotG = {cat(2,corrCondFrames{groups==max(groups)})};
    end
    cellfun(@(cm,cl) shadedErrorBar(1:size(cm,1),movmean(mean(cm,2,'omitnan'),3)',...
        movmean(std(cm,0,2,'omitnan'),3)./sqrt(length(dates)),'lineprops',{'Color',cl,'LineWidth',2}),...
        plotG,colorsG(3*(n-1) + 1 + double(n==max(groups)+1) : ...
        3*(n-1)+1 + sum(groups==min(n,max(groups)))-1 - double(n==max(groups)+1)));

    if(n==max(groups)+1)
        legend(flipud(ga.Children(isgraphics(ga.Children,'line'))), "[Rest][Movment]");
    else
        legend(flipud(ga.Children(isgraphics(ga.Children,'line'))),...
            condComNames(groups==min(n,max(groups)),1) + condComNames(groups==min(n,max(groups)),2),'Location','southwest');
    end
    ylim([-1 1]);
    saveas(gcf,['S:\Lab\ngc14\Working\',monkey,'\XCorr\',num2str(n)],'epsc');
    saveas(gcf,['S:\Lab\ngc14\Working\',monkey,'\XCorr\',num2str(n)],'png');
end