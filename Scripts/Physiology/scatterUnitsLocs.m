typeNames = ["Reach","Grasp","Both","All"];
conds = ["E","L","P"];
sOrder = [3 1 2];
pairs = combnk(1:length(sOrder),2);
load("S:\Lab\ngc14\Working\Both\Full_Baseline\FRTable_STransformed.mat");
repNames = ["Arm","Trunk","Hand"];
savePath = "S:\Lab\ngc14\Working\Both\Full_Baseline\SI\XY_Histograms\";
dir = ["ML","RC"];
cl = validatecolor(["#A2142F","#0072BD","#EDB120","#77AC30"],'multiple');
whiskerP = 10;
%%
close all;
monkeys = 1;
mStep = 55.65;
minBinPerc = 5;
newtPhys = tPhys;
sig = cell(length(dir),1);
allTable = table();
for c = 1:length(conds)
    %for r = 1:length(repNames)
    ec = {};
    ect = {};
    for d = 1:length(dir)
        sax = [];
        figure('Units','normalized','Position',[0 0 1 1]);
        figure('Units','normalized','Position',[0 0 1 1]);
        colormap(cell2mat(arrayfun(@(m) MotorMapping.repColors.(m),repNames,'UniformOutput',false)'));
        for m = 1:monkeys
            if(m==1)
                monkey= "Gilligan";
            else
                monkey = "Skipper";
            end
            pInd = [1:monkeys:length(typeNames)*monkeys]+(m-1);
            if(d==1)
                cBin = tPhys{:,'Y'};
            else
                cBin = tPhys{:,'X'};
            end
            tInd = tPhys.("TaskUnits_"+conds(c))>0;
            edgs = min(cBin(tInd)):mStep:(max(cBin(tInd))+(mStep-rem(max(cBin(tInd)),mStep)))+abs(rem(min(cBin(tInd)),mStep));
            [tvalue,~,~] = histcounts(cBin,edgs,'Normalization','count');
            totalUnits = sum(tvalue);
            tvalue = 100*(tvalue./totalUnits);
            startInd = find([-Inf tvalue]>minBinPerc,1,'first')-1;
            lastInd = find([-Inf tvalue]>minBinPerc,1,'last');
            tvalue(startInd) = sum(tvalue(1:startInd));
            tvalue(lastInd-1) = sum(tvalue(lastInd-1:end));
            tvalue(lastInd:end) = [];
            tvalue(1:startInd-1) = [];
            edgs(lastInd:end-1) = [];
            edgs(2:startInd) = [];
            bars = string(arrayfun(@num2str, (0:length(edgs)-1).*ImagingParameters.px2mm*mStep, 'UniformOutput',false));
            sig{d}(c,:,:) = zeros(1,length(typeNames),length(bars)-1);
            for s = 1:length(typeNames)
                if(s==length(typeNames))
                    sInd = tInd;
                else
                    sInd = tInd & tPhys{:,("unitType_E")}==s;
                    % ~isnan(tPhys.("rgSI_"+conds(c))) & [tPhys{:,"Somatotopy"}]==repNames(r) & strcmp(string(tPhys.Monkey),monkey)
                    % newtPhys = addvars(newtPhys,NaN(height(newtPhys),1),'NewVariableNames',svar,'After',width(newtPhys));
                end
                [hvalues,~,p] = histcounts(cBin,edgs,'Normalization','count');
                hvalues = hvalues./(sum(hvalues)); % (totalUnits*tvalue./100); %
                figure(1);
                sax(m,s) = subplot(monkeys,length(typeNames),pInd(s));
                title(sax(m,s),typeNames(s));
                hold on;
                if(0)
                    for ac = 1:(lastInd-startInd)
                        barReps = string(tPhys{sInd(p==unBars(ac)),'Somatotopy'});
                        [~,sortOrder] = ismember(barReps,repNames);
                        [~,sortOrder] = sort(sortOrder);
                        barReps = string(barReps(sortOrder));
                        tPoints = cumsum(arrayfun(@(m) sum(strcmp(barReps,m))/length(barReps), repNames));
                        patch([ac-1; ac; ac; ac-1], [0 0 hvalues(ac) hvalues(ac)], 'w','FaceColor',cl(s,:));
                        % hfl=patch([ac-1, ac-1, ac-1; ac, ac, ac; ac, ac, ac; ac-1, ac-1, ac-1],...
                        %     [0, hvalues(ac)*tPoints(1), hvalues(ac)*tPoints(2);0, hvalues(ac)*tPoints(1), hvalues(ac)*tPoints(2);...
                        %     hvalues(ac)*tPoints(1), hvalues(ac)*tPoints(2), hvalues(ac)*tPoints(3);...
                        %     hvalues(ac)*tPoints(1), hvalues(ac)*tPoints(2), hvalues(ac)*tPoints(3)],[1;2;3]);
                    end
                else
                    csBin = tPhys{sInd,"rgSI_"+conds(c)};
                    sig{d}(c,s,:) = accumarray(p(sInd),csBin,[max(p(sInd)),1],@ttest, 0,1);
                    if(s==4)
                        disp('');
                    end
                    boxplot(csBin,p(sInd),'Whisker',median(arrayfun(@(a) diff(prctile(csBin(p(sInd)==a),[25 75]+([1 -1].*whiskerP/2)))...
                        /iqr(csBin(p(sInd)==a)), 1:length(bars)-1)),'Notch','on','Widths',0.8,'Symbol','');
                    f = findobj(get(sax(m,s),'Children'),{'tag', 'Box', '-or','tag','Median'});
                    [f.Color] = deal(cl(s,:));
                    [f.LineWidth] = deal(2);
                    f = findobj(get(gca,'Children'),{'tag', 'Lower Whisker', '-or' ,'tag','Upper Whisker'});
                    [f.LineStyle] = deal('-');
                    [f.LineWidth] = deal(1.5);
                    %plot([0,max(pa)+1],[0 0],'LineWidth',1,'Color','k');
                end
                xticks([1:length(bars)]);
                xticklabels(strsplit(sprintf('%.2f ',bars(2:end))));
                xtickangle(25)
                if(d==1)
                    xlabel("Caudal/Rostral (microns)");
                else
                    xlabel("Medial/Lateral (microns)");
                end
                ylabel(sprintf("SI\n Grasp \t \t \t Reach"));
                yl = get(gca,'YLim');
                ylim([-1 1])
                plot([0, length(bars)],[0 0],'k');
                ec{s} = cBin(sInd);%.*sum(cell2mat(arrayfun(@(pp)double(p==pp).*(sum(p==pp)/length(p))*(sum(pt==pp)/length(pt)),unique(p),'UniformOutput',false)'),2);
                ect{s} = csBin;
            end
            figure(2);
            for si = 1:length(pairs)
                subplot(3,1,si);
                hold on;
                f1 = histcounts(categorical(round(ec{pairs(si,1)})),categorical(unique(round(ec{end}))));
                f2 = histcounts(categorical(round(ec{pairs(si,2)})),categorical(unique(round(ec{end}))));
                %f1=histcounts(ec{pairs(si,1)}.*ect{pairs(si,1)},length(unique(ec{end})));
                %f2=histcounts(ec{pairs(si,2)}.*ect{pairs(si,2)},length(unique(ec{end})));
                [e1,x1,u1,d1] = ecdf(unique(round(ec{end})),'Frequency',f1);
                pc = plot(x1,e1,'-',x1,u1,'--',x1,d1,'--');
                [pc.Color] = deal(cl(pairs(si,1),:));
                [pc.LineWidth] = deal(2);
                [e2,x2,u2,d2] = ecdf(unique(round(ec{end})),'Frequency',f2);
                pc = plot(x2,e2,'-',x2,u2,'--',x2,d2,'--');
                [pc.Color] = deal(cl(pairs(si,end),:));
                [pc.LineWidth] = deal(2);
                [~,pvalue(si),ks(si)] = kstest2(ec{pairs(si,1)},ec{pairs(si,2)});
                hs(si) = length(ec{pairs(si,1)}) + length(ec{pairs(si,2)});
                scatter([x1(find(e1>=.5,1,'first')), x2(find(e2>=.5,1,'first'))],[0 0], 120,[cl(pairs(si,1),:);cl(pairs(si,end),:)],'|','LineWidth',5)
                title(typeNames(pairs(si,1))+"/"+typeNames(pairs(si,2)));
                text(sum(get(gca,'XLim'))/2,0.8,"p = "+num2str(pvalue(si),'%.4f'),'Units','data');
                set(gca,'XTickLabels',str2double(get(gca,'XTickLabel')).*ImagingParameters.px2mm)
            end
            allTable(end+1:end+3,:) = cell2table([repmat({conds(c)},length(pvalue),1),num2cell(pairs,2),...
                repmat({dir(d)},length(pvalue),1),num2cell(pvalue)',num2cell(ks)'],...
                'VariableNames',{'cond','pair','axis','pvalue','ksvalue'});
            linkaxes([sax(m,1:end-1)]);
        end

        saveFigures(figure(1),savePath,"BothSites"+conds(c)+"_"+dir(d),[]);
        saveFigures(figure(2),savePath+"CDF\","BothSites"+conds(c)+"_"+dir(d),[]);
        close all;
    end
    %    end
end
writetable(allTable,savePath+"CDF\stats");