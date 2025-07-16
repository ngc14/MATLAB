monkey = 'Skipper';
strcmp(getenv('COMPUTERNAME'),'OMAR-ANALYSIS')
    conds = {'[ExtraSmallSphere]','[LargeSphere]', '[Photocell]','[Rest]'};
if(strcmp(monkey,'Gilligan'))
    mm = MotorMapping(42);
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
elseif(strcmp(monkey,'Skipper'))
    mm = MotorMapping(56);
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
end
LPk = 5;
HPk = 250;
refMask = imread(['S:\Lab\',monkey,'\Mapping\clean_mask_filled'],'bmp')<180;
refMask = refMask(:,:,1);
frames = containers.Map(conds, repmat({num2cell(1:70)},1,length(conds)));

ims = cell(1,length(conds));
dateLabels = {};
%%
for cn = 1:length(conds)
    currCond = conds{cn};
    allFrames = frames(currCond);
    imsF = {};
    dateTrials = {};
    [Otrans,spacing,allTrials] = deal(cell(1,length(dates)));
    for f = 1:length(allFrames)
        %%
        lf = allFrames{f};
        fi = clock;
        parfor (d = 1:length(dates))
            rPath = ['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},...
                '\Imaging\run0',num2str(runs{d}),'\Results\Blockview\',currCond,...
                '\FF_1\LP',num2str(LPk),'_HP',num2str(HPk),'\Align_Offset\noClip\'];
            if(exist([rPath, 'NaN_filtered_HP',num2str(HPk),'_LP',num2str(LPk),'.h5'],'file'))
                allTrials{d} = (squeeze(mean(h5read([rPath,'NaN_filtered_HP',...
                    num2str(HPk),'_LP',num2str(LPk),'.h5'],'/allTrials',...
                    [1 1 lf(1) 1],[size(refMask,1),size(refMask,2),...
                    length(lf),Inf]),3,'omitnan')));%,3,'omitnan');
                if(f==1)
                    tform = matfile(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},...
                        '\Imaging\run0',num2str(runs{d}),'\tform_nonrigid.mat']);
                    Otrans{d} = tform.O_trans;
                    spacing{d} = tform.Spacing;
                    tform = [];
                end
                for t = 1:size(allTrials{d},3)
                allTrials{d}(:,:,t) = bspline_transform(Otrans{d},allTrials{d}(:,:,t),spacing{d},3);
                end
            end
        end

        if(~exist(['S:\Lab\ngc14\Working\',monkey,'\',conds{cn},'\Frames\'],'dir'))
            mkdir(['S:\Lab\ngc14\Working\',monkey,'\',conds{cn},'\Frames\']);
        end
        imsF{f} = tall(mean(cat(3,allTrials{:}),3));
        %         imsF{f}(repmat(~allMMs{1}, [1 1 size(imsF{f},3)])) = NaN;
        %           condFrames{f} = mean(allTrials,3,'omitnan');
        %           stdCFrames{f} = std(allTrials,0,3,'omitnan');
        %     end
        %     save(['S:\Lab\ngc14\Working\',monkey,'\',conds{cn},'\MeanFrames.mat'],'condFrames')
        %     save(['S:\Lab\ngc14\Working\',monkey,'\',conds{cn},'\STDFrames.mat'],'stdCFrames');
        disp(strjoin(["Frame ",num2str(f)," done. (", num2str(etime(clock,fi),3)," s)"],''));
    end
    ims{cn} = imsF;%cell2mat(reshape(cellfun(@gather, imsF, 'UniformOutput',false), 1, 1, []))
    imsF = cell2mat(reshape(cellfun(@gather, ims{4}, 'UniformOutput',false), 1, 1, []));
    save(['S:\Lab\ngc14\Working\', monkey,'\', currCond,'\TrialAVG'], 'imsF','-v7.3')
end
%%
% saveas(gcf,strjoin(["S:\Lab\ngc14\Working\",monkey,"\Corr\Forelimb\All_Corrs_Clipped.png"],''));