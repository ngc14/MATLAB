function [imOut] = imageFilter_LPHP_ngc14_NaN(imageBig,LPk,HPk,mask)
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
if HPk > 0 && ~isempty(imageBig)
    HPk = round(HPk/4) + double(rem(round(HPk/4),2)==0);
    imageHP = imageBigLP;
    imageHP(~mask) = NaN;

    if(numel(size(imageBig))==3)
        imageHP = imresize3(imageHP,[fix(size(imageHP,[1 2]).*.25),size(imageHP,3)],'Method','linear');
        maskSmall = imresize3(mask,[fix(size(imageHP,[1 2]).*.25),size(imageHP,3)],'Method','linear');

        %imageF = medfilt3(imageHP,[HPk,HPk,1],'symmetric');
        imageP = NaN(size(imageHP,1)+HPk, size(imageHP,2)+HPk,size(imageHP,3));
        imageP(round(HPk/2):end-round(HPk/2),round(HPk/2):end-round(HPk/2),:) = imageHP;
        imageP = parallel.pool.Constant(imageP);
        xSz = size(imageHP,1);
        ySz = size(imageHP,2);
        [xn,yn] = ind2sub(size(imageHP,[1 2]),1:(xSz*ySz));
        xyWin = reshape(arrayfun(@(xs,xe,y,ye) {xs:xe, y:ye}, ...
            xn,xn+HPk,yn,yn+HPk,'UniformOutput',false),size(imageHP,[1 2]));
        parfor x = 1:xSz
            for y = 1:ySz
                imageF(x,y,:) = median(imageP.Value(xyWin{x,y}{1},xyWin{x,y}{2},:),[1 2],'omitnan');
            end
        end

        imageF(isnan(imageF)) = median(imageF(maskSmall),'all', 'omitnan');
        imageBigHP = imresize3(gather(imageF),size(imageBig),'Method','linear');
    else
        imageHP = imresize(imageHP,.25,'Method','bilinear');
        imageF = medfilt2(imageHP,[HPk HPk]);
        imageBigHP = imresize(imageF,size(imageBig),'Method','bilinear');
    end
    imOut = imageBigLP - imageBigHP;
else
    imOut = imageBigLP;
end
end