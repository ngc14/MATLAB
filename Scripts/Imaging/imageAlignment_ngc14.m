 function [] = imageAlignment_ngc14(varargin)
%% user input

username = 'ngc14';

animal = 'Skipper';
date = '11_06_2020';
run = 'run00';
stimSelect = [0 2 4];

if ~isempty(varargin)
    % dates
    if ~isempty(varargin{1})
        date = varargin{1};
    end
    
    % runs
    if ~isempty(varargin{2})
        run = ['run0',num2str(varargin{2})];
    end
    
    % stims
    if ~isempty(varargin{3})
        stimSelect = varargin{3};
    end
end


alignType = 'similarity';   % https://www.mathworks.com/help/images/ref/fitgeotrans.html
% https://www.mathworks.com/examples/computer-vision/mw/vision-ex57604986-image-registration-using-multiple-features

threshold = 10;
%% variable initialization, find block data
parent = ['\\pitt\sni\gharbawie\Lab\',username,'\Scripts\'];
path = ['S:\Lab\',animal,'\All Data\',animal,'_',date, '\Imaging\',run];

cd(path);
blocks = dir('*.BLK');
[~,indx] = natsort({blocks.name});
blocks = blocks(indx);

anapar = OIHeadRead((blocks(1).name), 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
NBlocks = length(blocks);
conds = [];
for n =1:size(blocks,1)
    anapar = OIHeadRead((blocks(n).name), 'v');
    conds = [conds, anapar.Cond];
end
numConds = unique(conds);
timeStart = clock;
dispstat('Loading & Aligning Data... 0.00% Done','init');
% addAttachedFiles(gcp,'C:\Program Files\MATLAB\R2018a\toolbox\vision\vision\detectBRISKFeatures.m');
% addAttachedFiles(gcp,'C:\Program Files\MATLAB\R2018a\toolbox\vision\vision\detectSURFFeatures.m');
% addAttachedFiles(gcp,'C:\Program Files\MATLAB\R2018a\toolbox\vision\vision\extractFeatures.m');
% addAttachedFiles(gcp,'C:\Program Files\MATLAB\R2018a\toolbox\vision\vision\matchFeatures.m');
% addAttachedFiles(gcp,'C:\Program Files\MATLAB\R2018a\toolbox\vision\vision\estimateGeometricTransform.m');
%addAttachedFiles(gcp,'\\pitt\sni\gharbawie\Lab\ngc14\Scripts\Imaging\');

% align to first frame of first block
frames = OIReadStim(blocks(round(length(blocks)/2)).name, 1, 'v');
original = frames(:,:,1);
original = original - mean(original(:)); % mean-subtract and clip image to increase contrast
[original, ~, clipLow, clipHigh] = OIClipH2(original, 1, 1.5, ones(H,W));
frame_original = original;

% detect features in original image
ptsOriginalBRISK  = detectBRISKFeatures(original, 'MinContrast', 0.9,'MinQuality',0.9);
ptsOriginalSURF  = detectSURFFeatures(original);
[featuresOriginalFREAK,  validPtsOriginalBRISK]  = extractFeatures(original,  ptsOriginalBRISK);
[featuresOriginalSURF,  validPtsOriginalSURF]  = extractFeatures(original,  ptsOriginalSURF);
outputView = imref2d(size(original));
%% main alignment loop
for c = 1:length(stimSelect)
    condBlocks = blocks(find(conds==numConds(c)));
    distCond = cell(size(condBlocks,1), NFrames);
    displacement_xCond =cell(size(condBlocks,1), NFrames);
    displacement_yCond = cell(size(condBlocks,1), NFrames);
    pts_originalCond = cell(size(condBlocks,1), NFrames);
    pts_distortedCond = cell(size(condBlocks,1), NFrames);
    for n = 1:size(condBlocks,1)
        loopStart = clock;
        frames = OIReadStim(condBlocks(n).name, 0, 'v');
        
        % align all frames
        parfor f = 1:NFrames
            % prepare images for alignment
            distorted = frames(:,:,f);
            distorted = distorted - mean(distorted(:)); % mean-subtract and clip image to increase contrast
            [distorted, ~, ~, ~] = OIClipH2(distorted, 9, [clipLow clipHigh], ones(H,W));
            
            % use BRISK and SURF for image alignment
            ptsDistortedBRISK = detectBRISKFeatures(distorted, 'MinContrast', 0.9,'MinQuality',0.9);
            ptsDistortedSURF = detectSURFFeatures(distorted);
            [featuresDistortedFREAK, validPtsDistortedBRISK] = extractFeatures(distorted, ptsDistortedBRISK);
            [featuresDistortedSURF, validPtsDistortedSURF] = extractFeatures(distorted, ptsDistortedSURF);
            
            try
                indexPairsBRISK = matchFeatures(featuresOriginalFREAK, featuresDistortedFREAK, 'MatchThreshold', 40, 'MaxRatio', 0.8);
                indexPairsSURF = matchFeatures(featuresOriginalSURF, featuresDistortedSURF);
                matchedDistortedBRISK = validPtsDistortedBRISK(indexPairsBRISK(:,2));
                matchedDistortedSURF = validPtsDistortedSURF(indexPairsSURF(:,2));
                matchedDistortedXY = [matchedDistortedSURF.Location; matchedDistortedBRISK.Location];
                matchedOriginalBRISK  = validPtsOriginalBRISK(indexPairsBRISK(:,1));
                matchedOriginalSURF  = validPtsOriginalSURF(indexPairsSURF(:,1));
                matchedOriginalXY  = [matchedOriginalSURF.Location; matchedOriginalBRISK.Location];
                
                
                % initial alignment
                [tformTotal,inlierDistortedXY,inlierOriginalXY] = ...
                    estimateGeometricTransform(matchedDistortedXY,matchedOriginalXY,alignType,...
                    'maxNumTrials',10000,'confidence',99,'MaxDistance',5);
                
                % get rid of outlier point pairs
                distTemp = sqrt((inlierOriginalXY(:,1)-inlierDistortedXY(:,1)).^2 + ...
                    (inlierOriginalXY(:,2)-inlierDistortedXY(:,2)).^2);
                distThreshHigh = median(distTemp) + 3*std(distTemp);
                distThreshLow = median(distTemp) - 3*std(distTemp);
                
                if std(distTemp) > 0
                    inlierDistortedXY = inlierDistortedXY((distTemp > distThreshLow) & (distTemp < distThreshHigh),:);
                    inlierOriginalXY = inlierOriginalXY((distTemp > distThreshLow) & (distTemp < distThreshHigh),:);
                    
                    % re-do alignment without outlier points
                    [tformTotal,inlierDistortedXY,inlierOriginalXY] = ...
                        estimateGeometricTransform(inlierDistortedXY,inlierOriginalXY,alignType,...
                        'maxNumTrials',10000,'confidence',99,'MaxDistance',5);
                    
                    distCond{n,f} = sqrt((inlierOriginalXY(:,1)-inlierDistortedXY(:,1)).^2 + ...
                        (inlierOriginalXY(:,2)-inlierDistortedXY(:,2)).^2);
                else
                    distCond{n,f} = distTemp;
                end
                
                % x/y displacement between matched point pairs
                displacement_xCond{n,f} = inlierDistortedXY(:,1)-inlierOriginalXY(:,1);
                displacement_yCond{n,f} = inlierDistortedXY(:,2)-inlierOriginalXY(:,2);
                
                %frame_distorted{n,f} = distorted;
                pts_originalCond{n,f} = inlierOriginalXY;
                pts_distortedCond{n,f} = inlierDistortedXY;
                %align_XYZRS(:,:,n,f) = tformTotal.T;
            catch
                disp('Could not match features');
                distCond{n,f} = NaN;
                % x/y displacement between matched point pairs
                displacement_xCond{n,f} = NaN;
                displacement_yCond{n,f} = NaN;
                
                %frame_distorted{n,f} = distorted;
                pts_originalCond{n,f} = NaN;
                pts_distortedCond{n,f} = NaN;
            end
        end
        
        % display function progress
        percDone = (c-1)/length(numConds) + (1/length(numConds) * n/size(condBlocks,1));
        percDone = round(percDone*10000)/100;
        elapsedTime = etime(clock,timeStart);
        dispstat(['Loading & Aligning Data... ',num2str(percDone),'% Done.']);
    end
    dist{c} = distCond;
    displacement_x{c} = displacement_xCond;
    displacement_y{c} = displacement_yCond;
    pts_original{c} = pts_originalCond;
    pts_distorted{c} = pts_distortedCond;
end

totalTime = round(etime(clock,timeStart)*100/60)/100;
dispstat(['Loading & Aligning Data... 100.00% Done. Total time: ',num2str(totalTime),' min.']);

%cd(parent);

%% movie 2
if 1==0
    figure;
    for n = 1:NBlocks
        frames = OIReadStim(blocks(n).name, stimSelect, 'v');
        for f = 1:NFrames
            
            distorted = frames(:,:,f);
            distorted = distorted - mean(distorted(:)); % mean-subtract and clip image to increase contrast
            [distorted, ~, ~, ~] = OIClipH(distorted, 9, [clipLow clipHigh], ones(H,W));
            imagesc(distorted);
            axis image off; colormap gray;
            title(['Block ',num2str(n),' Frame ',num2str(f)]);
            
            pause(0.1);
        end
    end
end


%% distance box plot
f1 = figure;
f2 = figure;
for c = 1:length(numConds)
    set(0,'CurrentFigure', f1);
    subplot(3,length(numConds),c);
    if(c==1)
        ylabel('Median Dist Offset (px)');
        xlabel('Block #');
    end
    title(numConds(c))
    medVals = cellfun(@nanmedian, dist{c})';
    stdVals = cellfun(@nanstd, dist{c})';
    shadedErrorBar(1:1/NFrames:size(dist{c},1)+(1-1/NFrames),...
        medVals(:),stdVals(:),'lineprops',{'-','linewidth',0.5,'color','r'})
    
    set(0,'CurrentFigure', f2);
    subplot(3,length(numConds),c);
    if(c==1)
        xlabel('Frame #');
        ylabel('Median Dist Offset (px)');
    end
    title(numConds(c))
    hold on;
    for n = 1:size(dist{c},1)
        plot(1:NFrames, cellfun(@nanmedian, dist{c}(n,:)), 'r-');
    end
    plot(1:NFrames, median(cellfun(@median, dist{c})), 'k-', 'LineWidth', 2);
    hold off;
end
for c = 1:length(numConds)
    set(0,'CurrentFigure', f1);
    subplot(3,length(numConds),length(numConds)+c);
    if(c==1)
        xlabel('Block #');
        ylabel('Median X Offset (px)');
    end
    title(numConds(c))
    medVals = cellfun(@nanmedian, displacement_x{c})';
    stdVals = cellfun(@nanstd, displacement_x{c})';
    shadedErrorBar(1:1/NFrames:size(displacement_x{c},1)+(1-1/NFrames),...
        medVals(:),stdVals(:),'lineprops',{'-','linewidth',0.5,'color','g'});
    
    set(0,'CurrentFigure', f2);
    subplot(3,length(numConds),length(numConds)+c);
    if(c==1)
        xlabel('Frame #');
        ylabel('Median X Offset (px)');
    end
    title(numConds(c))
    hold on;
    for n = 1:size(dist{c},1)
        plot(1:NFrames, cellfun(@nanmedian, displacement_x{c}(n,:)), 'g-');
    end
    plot(1:NFrames, median(cellfun(@median, displacement_x{c})), 'k-', 'LineWidth', 2);
end


for c = 1:length(numConds)
    set(0,'CurrentFigure', f1);
    subplot(3,length(numConds),2*length(numConds)+c);
    if(c==1)
        xlabel('Block #');
        ylabel('Median Y Offset (px)');
    end
    title(numConds(c))
    medVals = cellfun(@nanmedian, displacement_y{c})';
    stdVals = cellfun(@nanstd, displacement_y{c})';
    shadedErrorBar(1:1/NFrames:size(displacement_y{c},1)+(1-1/NFrames),...
        medVals(:),stdVals(:),'lineprops',{'-','linewidth',0.5,'color','b'});
    
    set(0,'CurrentFigure', f2);
    subplot(3,length(numConds),2*length(numConds)+c);
    if(c==1)
        xlabel('Frame #');
        ylabel('Median Y Offset (px)');
    end
    title(numConds(c))
    hold on;
    for n = 1:size(dist{c},1)
        plot(1:NFrames, cellfun(@nanmedian, displacement_y{c}(n,:)), 'b-');
    end
    plot(1:NFrames, median(cellfun(@median, displacement_y{c})), 'k-', 'LineWidth', 2);
end
disp_x_med = cellfun(@(a) cellfun(@nanmedian, a, 'UniformOutput', false), displacement_x, 'UniformOutput', false);
disp_y_med = cellfun(@(a) cellfun(@nanmedian, a, 'UniformOutput', false), displacement_y, 'UniformOutput', false);
for c = 1:length(numConds)
    [row,col] =find(cellfun(@(a) median(a)>threshold, dist{c}));
    row = row(~(ismember(col,11:20)));
    badAlign{c} = unique(row);
end
%%
saveas(f1, 'Blocks.fig')
saveas(f2, 'Trials.fig')
save('alignment_offset.mat', 'disp_x_med', 'disp_y_med');
save('offsets.mat', 'displacement_x', 'displacement_y', 'dist');
save('bad_alignments.mat', 'badAlign');

% end