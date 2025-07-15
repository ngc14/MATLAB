
%% user variables
regressSkull = false;            % regress skull?
    numSkullMasks = 4;          % how many skull masks?
regressVessel = true;           % regress vessel?
    numVesselMasks = 2;         % how many vessel masks?
regressMean = false;            % regress average image signal?
regressorSmoothKernel = 5;      % kernel for moving median smooth on regressor signals
LOAD_REGRESSION_MASKS = false;  % load pre-defined vessel masks?


%% load or select masks
if regressSkull && ~LOAD_REGRESSION_MASKS   % select skull masks
    figure('position',[300 100 1200 800]);
    
    show(greenGraySmall);
    title('Select Skull Mask');

    for i = 1:numSkullMasks
        skullMaskRR(:,:,i) = roipoly;
        hold on;
        imagesc(ones(size(skullMaskRR,1),size(skullMaskRR,2),3),'alphadata',skullMaskRR(:,:,i)*0.5);
    end
    
    show(im_super(greenGraySmall,sum(skullMaskRR,3),0.3));
    title('Skull Mask');
    pause(0.1);
    save(['regMaskSkull_SS',num2str(SPATIAL_BIN),'.mat'],'skullMaskRR');
end

if regressVessel && ~LOAD_REGRESSION_MASKS  % select vessel masks
    
    title('Select Vessel Mask');

    for i = 1:numVesselMasks
        vesselMaskRR(:,:,i) = roipoly;
        hold on;
        imagesc(ones(size(vesselMaskRR,1),size(vesselMaskRR,2),3),'alphadata',vesselMaskRR(:,:,i)*0.5);
    end
    
%     show(im_super(greenGraySmall,(sum(vesselMaskRR,3)+sum(skullMaskRR,3))/2,0.3));
    title('Regression Masks');
    pause(0.1);
    save(['regMaskVessel_SS',num2str(SPATIAL_BIN),'.mat'],'vesselMaskRR');
end


if LOAD_REGRESSION_MASKS % load predefined masks?
    load(['regMaskSkull_SS',num2str(SPATIAL_BIN),'.mat']);
    load(['regMaskVessel_SS',num2str(SPATIAL_BIN),'.mat']);
end



%% perform regression after loading data
if regressSkull || regressVessel || regressMean

    dispstat('','init');
    dispstat('PERFORMING REGRESSION...');
    
    for i = 1:size(frames_2Hz,3) % get average regression mask values at each time step
        tempFrame = frames_2Hz(:,:,i);
        
        if regressSkull
            for j = 1:numSkullMasks 
                skullSig(i,j) = mean(tempFrame(skullMaskRR(:,:,j)));
            end
        end

        if regressVessel
            for j = 1:numVesselMasks
                vesselSig(i,j) = mean(tempFrame(vesselMaskRR(:,:,j)));
            end
        end
    end

    % smooth out regression signals a little bit
    if regressSkull, skullSig = movmedian(skullSig,regressorSmoothKernel); end
    if regressVessel, vesselSig = movmedian(vesselSig,regressorSmoothKernel); end
    
    % reshape 3d image matricies into 2d ones
    frames_2Hz_reshape = reshape(frames_2Hz,[size(frames_2Hz,1)*size(frames_2Hz,2) size(frames_2Hz,3)])';
    
    % construct regression matrix X
    X = [ones(size(frames_2Hz_reshape,1),1)];
    
    if regressSkull
        for j = 1:numSkullMasks
            X = [X, skullSig(:,j)];
        end
    end
    
    if regressVessel
        for j = 1:numVesselMasks
            X = [X, vesselSig(:,j)];
        end
    end
    
    if regressMean
        X = [X, mean(frames_2Hz_reshape,2)];
    end
    
    % perform regression
    for i = 1:size(frames_2Hz_reshape,2)
        b(i,:) = regress(frames_2Hz_reshape(:,i),X);
        frames_2Hz_reshape_reg(:,i) = frames_2Hz_reshape(:,i)-X*b(i,:)';
    end
    
    % reshape back into 3d image
    frames_2Hz_reg = reshape(frames_2Hz_reshape_reg',[size(frames_2Hz,1) size(frames_2Hz,2) size(frames_2Hz,3)]);
    for i=1:size(b,2)
        br{i} = reshape(b(:,i),[size(frames_2Hz,1) size(frames_2Hz,2)]);
    end
%     frames_2Hz = frames_2Hz_reg;

    clear frames_2Hz_reshape*; 
    dispstat('PERFORMING REGRESSION... 100% Done.');
end
