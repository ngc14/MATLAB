function [imOut] = imageFilter_LPHP_ngc14_NaN_GPU(imageBig,LPk,HPk,mask)
if isempty(mask)
    mask = true(size(imageBig));
end
%% low pass
if LPk>0 && ~isempty(imageBig)
    if(numel(size(imageBig))==3)
        mask2 = imgaussfilt3(double(mask),round(LPk/2),'FilterSize', [LPk, LPk, 1], 'Padding', 0,'FilterDomain', 'spatial');
        A_masked = imgaussfilt3(imageBig.*mask,round(LPk/2),'FilterSize', [LPk, LPk, 1], 'Padding', 0,'FilterDomain', 'spatial');
    else
        mask2 = imgaussfilt(double(mask),round(LPk/2),'FilterSize', [LPk, LPk], 'Padding', 0,'FilterDomain', 'spatial');
        A_masked = imgaussfilt(imageBig.*mask,round(LPk/2),'FilterSize', [LPk, LPk], 'Padding', 0,'FilterDomain', 'spatial');
    end
    mask3 = mask2+double(mask2 == 0);
    A_masked = A_masked./mask3;
    imageBigLP = A_masked.*double(mask2 ~= 0) + imageBig.*double(mask2 == 0);
    imageBigLP(mask==0) = imageBig(mask==0);
else
    imageBigLP = imageBig;
end
%% high pass
HPk = round(HPk/4) + double(rem(round(HPk/4),2)==0);
imageHP = imageBigLP;
imageHP(~mask) = NaN;
if(numel(size(imageBig))==3)
    imageHP = imresize3(imageHP,[fix(size(imageHP,[1 2]).*.25),size(imageHP,3)],'Method','linear');
    maskSmall = imresize3(mask,size(imageHP),'Method','linear');
    imageF = cell(size(imageHP,[1 2]));
else
    imageHP = imresize(imageHP,.25,'Method','bilinear');
    imageF = medfilt2(gpuArray(imageHP),[HPk HPk]);
    imageBigHP = imresize(imageF,size(imageBig),'Method','bilinear');
end
imageP = NaN(size(imageHP,1)+HPk, size(imageHP,2)+HPk,size(imageHP,3));
imageP(HPk/2+1:end-HPk/2,HPk/2+1:end-HPk/2,:) = imageHP;
[rn,cn] = ind2sub(size(imageHP,[1 2]),1:prod(size(imageHP,[1 2])));
imageP = gpuArray(imageP(:));

    function returnedPx = medFiltK(r,c,HPK)
        returnedPx = squeeze(median(imageP(r:r+HPK,c:c+HPK,:),[1 2],'omitnan'));
    end
for li = 1:length(rn)*double(HPk>0)
    imageF{li} = arrayfun(@medFiltK,rn(li),cn(li),HPk,'UniformOutput',false);
    disp(li);
end

if(numel(size(imageBig))==3)
    imageF = cell2mat(gather(imageF));
    imageF(isnan(imageF)) = median(imageF(maskSmall),'all','omitnan');
    imageBigHP = imresize3(gather(imageF),size(imageBig),'Method','linear');
end

imOut = imageBigLP - imageBigHP;
end