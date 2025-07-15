close all;
clear all;

monkey = 'Gilligan';                   % animal name
domain = 'M1';
condition = '[ExtraSmallSphere]-[BASELINE]';
pValue = 'p-0.0001';
frames = [51:56];
LPk = 15;
HPk = 550;
saveVar = 0;
numSessions = 0;

if(strcmp(monkey, 'Gilligan'))
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
elseif(strcmp(monkey,'Skipper'))
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0],[0],[0],[0]};
end

refMask = imread(['S:\Lab\',monkey,'\Mapping\clean_mask_filled'], 'bmp')>200;
if(strcmp(monkey,'Skipper'))
    refMask = refMask(:,2:end-1,:);
end

% for frameRanges = 1:1:70
%     frames = max(1,frameRanges-2):min(70,frameRanges+2);
    imT = zeros(size(refMask,1),size(refMask,2));
    for d = 1:length(dates)
        sessionPath = ['S:\Lab\',monkey,'\All Data\', monkey, '_',dates{d}, '\Imaging\'];
        % get largest run folder
        runFolder = [sessionPath,'run0', num2str(runs{d})];
        
        tTestPath = [runFolder, '\Results\Ttest_nsc15_V3\Frames',...
            num2str(frames(1)),'-',num2str(frames(end)),'\LP', num2str(LPk), 'HP',...
            num2str(HPk),'\Align_Offset\FF_1\Warped'];
        filesInPath = dir(tTestPath);
        if(any(cellfun(@(a) contains(a,condition), {filesInPath.name})))
            bmpIm = imread([tTestPath,'\NONRIGID_', condition, '_', pValue, '.bmp']);
            %imshow(bmpIm);
            if(~exist('imT','var'))
                imT = zeros(size(bmpIm));
            end
            imT = imT + bmpIm;
            numSessions = numSessions + 1;
        end
        %pause();
     end
    
    disp(max(max(imT)));
    figure();
    h = imshow(imT);
    name = [condition(isletter(condition) | condition=='-'), '_',pValue];
    %title([name,', ', num2str(numSessions),' sessions']);
    colorbar;
    set(h, 'alphadata', imT ~= 0);
    colormap jet;
    caxis([0 4]);
    colorbar off;
    hold on;
    h = imagesc(refMask);
    set(h, 'alphadata', refMask(:,:,1)==0);
    if(saveVar)
        save(['S:\Lab\', monkey,'\Mapping\tTests\', name, '_', num2str(frames),'.mat'],'imT');
    end
%end