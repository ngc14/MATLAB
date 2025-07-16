typeNames = ["Both","Reach", "Grasp"];
mk = ["o","o","square","x"];
conds = ["E","L","P"];
sOrder = [3 1 2];
savePath = "S:\Lab\ngc14\Working\Both\Baseline_FR\";
load(["X"+savePath{1}(2:end)+"FRTable.mat"]);
freps = string(fieldnames(MotorMapping.repColors));
freps = freps(contains(freps,string(tPhys.Somatotopy)));

mtform = matfile("S:\Lab\ngc14\Working\S2G_Reference_tform.mat");
gi = im2gray(imread(["S:\Lab\Gilligan\Mapping\clean_reference_mask.png"],"png",'BackgroundColor',[1 1 1]));
sk = im2gray(imread(["S:\Lab\Skipper\Mapping\clean_reference_mask.png"],"png",'BackgroundColor',[1 1 1]));
skc = imcrop(imwarp(sk,mtform.tform),[695 305 450 650]);
gkc = imcrop(gi,[1 1 450 650]);
%%
close all;
for c = 1:length(conds)
    f1 = figure('Units','normalized','Position',[0 0 1 1 ]);
    sax = [];
    for m = 1:1
        mPerc = 1.5;
        mSteps = 3;
        mxLim = [0 365];
        myLim = [0 600];
        if(m==1)
            monkey= "Gilligan";
            mSz = 85;
            boundaries = cell2mat(bwboundaries(im2gray(imread(["S:\Lab\"+monkey+"\Mapping\clean_reference_mask.png"],...
                "png",'BackgroundColor',[1 1 1])),"noholes"));
        else
            monkey = "Skipper";
            mSz = 125;
            boundaries = cell2mat(bwboundaries(imcrop(imwarp(im2gray(imread(["S:\Lab\"+monkey+"\Mapping\clean_reference_mask.png"],...
                "png",'BackgroundColor',[1 1 1])),mtform.tform),[695 305 450 650]),"noholes"));
        end
        mSz = 100; boundaries = cell2mat(bwboundaries(~(gkc<225 | skc<225),'noholes'));
        iminf = imfinfo(["S:\Lab\"+monkey+"\Mapping\clean_reference_mask.png"]);
        gf = ones(iminf.Width,iminf.Height);
        for a = 1:length(boundaries)
            gf(boundaries(a,1),boundaries(a,2)) = 0;
        end
        gf(1:10,:) = 1;
        gf(end-10:end,:) = 1;
        gf(:,1:1) = 1;
        gf(:,end-10:end) = 1;
        boundaries = bwboundaries(imfill(imdilate(~logical(gf),ones(10,10),'same'),'holes'),4,'TraceStyle','pixeledge');
        for s = 1:3
            tInd = tPhys.("TaskUnits_"+conds(c))>0 %& strcmp(string(tPhys.Monkey),monkey);
            %[vals,xEdges,yEdges,binX,binY] = histcounts2(tPhys{tInd,'Y'},tPhys{tInd,'X'},[6,5],'Normalization','percentage');
            binX = tPhys{tInd,'X'};
            binY = tPhys{tInd,'Y'};
            fInd = unique([binX,binY],'rows');
            cl = NaN(length(fInd),1);
            al = strings(length(fInd),1);
            for sp = 1:length(fInd)
                cl(sp) = (100*sum(tPhys{all(tPhys{:,{'X','Y'}}==fInd(sp,:),2) & tInd,("unitType_"+conds(c))}==...
                    sOrder(mod(s-1,3)+1)))/sum(tPhys{tInd,("unitType_"+conds(c))}==sOrder(mod(s-1,3)+1));
                currLabs = tPhys{all(tPhys{:,{'X','Y'}}==fInd(sp,:),2) & tInd,'Somatotopy'};
                if(~isempty(currLabs))
                    neum = dictionary(string(unique(currLabs)),arrayfun(@(l) sum(currLabs==l)/length(currLabs), unique(currLabs)));
                    ks = neum.keys('cell');
                    if(sum(neum.values==neum(mode(currLabs)))>1)
                        al(sp) = strjoin(string(ks(cellfun(@(k) neum(k) == neum(mode(currLabs)), ks))),",");
                    else
                        al(sp) = string(mode(currLabs));
                    end
                    % text(yBinCenters(yc),xBinCenters(xc),10,"["+strjoin(arrayfun(@(a,s) a{1}(1)+":"+ round(s)+"% ",...
                    %     string(unique(currLabs)), 100.* neum.values, 'UniformOutput',true))+"]",'HorizontalAlignment','center', ...
                    % 'VerticalAlignment','bottom','Rotation',45);
                end
            end
            %scatterbar3(reshape(c,size(vals)),reshape(r,size(vals)),vals,1,cl,al);
            cm = interp1(linspace(0,1,256),[interp1(linspace(0,1,2), ...
                [1 .85 .85; .93 .07 .07], linspace(0,1,255));.45 0 0], linspace(0,1,mSteps));
            set(0,'CurrentFigure', f1);
            sax(m,s) = subplot(1,3,(3*(m>1))+s);
            hold on;
            if(s==1 && m==1)
                for mm = 1:length(mk)-1
                    scatter(NaN,NaN,mSz,'k','filled',mk(mm));
                end
                legend(freps,'Location','southeast','AutoUpdate','off')
            end
            cellfun(@(b) plot(sax(m,s),b(:,2), b(:,1), 'Color', [0 .15 0], 'LineWidth',4), boundaries);
            clInd = discretize(cl,linspace(0,mPerc,mSteps));
            clInd(isnan(clInd)) = mSteps;
            clm = cm(clInd,:);
            clm(cl==0,:) = repmat([.70 .75 .75],sum(cl==0),1);
            colormap(cm);
            for ss = 1:length(fInd)
                sc = scatter(fInd(ss,1),fInd(ss,2),mSz,clm(ss,:),'o','filled');
                if(~contains(al(ss),","))
                    set(sc,'Marker',mk(strcmp(freps,al(ss))));
                    set(sc,'MarkerEdgeColor',[0.5 0.5 0.5]);
                else
                    set(sc,'MarkerEdgeColor',get(sc,'CData'));
                    set(sc,'MarkerFaceColor','none');
                    set(sc,'Marker',mk(end));
                end
                set(sc,'LineWidth',.5);
            end
            colorbar(sax(m,s),'eastoutside','Visible','on','Ticks',linspace(0,mPerc,mSteps+1));
            axis image
            title(sax(m,s),typeNames(mod(s-1,3)+1));
            set(gca,'YDir','reverse');
            set(gca,'XDir','normal');
            xlim(mxLim);
            ylim(myLim);
            xticklabels(20.*(get(gca,'XTick')));
            yticklabels(20.*(get(gca,'YTick')));
            %view(15,25);
            %xlb.Position = xlb.Position + [0 0 0];
            %ylb.Position = ylb.Position - [-.5 3 0];
            xlb = xlabel("Caudal/Rostral (microns)");%,'Rotation',-10;
            ylb = ylabel("Lateral/Medial (microns)");%,'Rotation',-110;
            clim([0 mPerc]);
        end
        linkaxes([sax(m,:)])
    end
    saveFigures(gcf,savePath+"Misc\Grids\","BothSites"+conds(c),[]);
end