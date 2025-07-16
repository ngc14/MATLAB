dirName = '\\130.49.229.252\gharbawie\Lab\Gilligan\Sorted EMG Data\Results\Models\';
models = dir(fullfile(dirName, '*.mat'));
modelNames = {'\bfBicep.mat', '\bfTricep.mat', '\bfDeltoid.mat', 'D4D5Extensor.mat', 'Extensor.mat', ...
    'ProximalRadialFlexor.mat', 'DistalFlexor.mat', 'ProximalUlnarFlexor'};
currModel = load([dirName, models(1).name]);
currModel = currModel.derivativesModel;
conds = fieldnames(currModel);
for c = 1:size(conds,1)
    condCorrMatrix = [];
    for m = 1:size(models,1)
        currModel = load([dirName, models(m).name]);
        currModel = currModel.derivativesModel.(conds{c});
        for i = 1:size(models,1)
            compModel = load([dirName, models(i).name]);
            compModel = compModel.derivativesModel.(conds{c});
            for s = 1:size(currModel,2)
                averageDerivModel = nanmean(cell2mat(currModel(:,s)));
                % photocell
                if(~isempty(averageDerivModel))
                    condCorrMatrix(m,i,s) = corr(averageDerivModel', nanmean(cell2mat(compModel(:,s)),1)');
                end
            end
        end
    end
    corrMatrix.(conds{c}) = condCorrMatrix;
end
%% Plot
n = 100;                %// number of colors

colors(1,:) = [1 1 1];   %// color first row - red
colors(2,:) = [0 0 1];   %// color 25th row - green

[X,Y] = meshgrid(1:3,1:n);  %// mesh of indices

cmap = interp2(X([1,n],:),Y([1,n],:),colors,X,Y); %// interpolate colormap

segs = {'Reach', 'Lift', 'Withdraw'};

for s = 1:3
    figure('Name', segs{s});
    for c = 1:size(conds,1)
        subplot(2,2,c);
        currCorrMat = corrMatrix.(conds{c});
        currCorrMat = currCorrMat(:,:,s);
        currCorrMat = max(currCorrMat,0);
        for r = 1:size(currCorrMat,1)
            for col = 1:size(currCorrMat,2)
                if(col>r)
                    currCorrMat(r,col) = NaN;
                end
            end
        end
        im = imagesc(currCorrMat);
        colormap(cmap) %// set color map
        caxis([0 1]);
        title(conds{c});
        ax = gca();
        ax.XTickLabel = modelNames;
        ax.YTickLabel = modelNames;
        ax.XTickLabelRotation = 30;
        line([3.5 3.5], [0 8.5], 'LineWidth', 2, 'Color', 'k');
        line([0 8.5], [3.5 3.5], 'LineWidth', 2, 'Color', 'k');
        set(im,'alphadata',~isnan(currCorrMat));
        set(ax,'color','black');
    end
end