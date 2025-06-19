date = '12_12_2018';
pVal = '0.001';
cond = 'LargeSphere';
pathName = ['\\pitt\sni\Gharbawie\Lab\Gilligan\All Data\Gilligan_',date,...
    '\Imaging\run00\Results\tTest\',cond,'\LP10_HP700\p-',pVal];
numCols = 10;

direc = dir(pathName);
direc = direc(cellfun(@(a) contains(a, '.bmp'), {direc.name}));
direc = direc(cellfun(@(a) ~contains(a, 'tTest'), {direc.name}));
natsort({direc.name});
[~,in] =natsort({direc.name});
direc = direc(in);
green = imread(['\\pitt\sni\Gharbawie\Lab\Gilligan\All Data\Gilligan_',date,...
    '\Imaging\run00\green00_edited.bmp']);
green = rgb2gray(green);
[H,W] = size(green);
NFrames = size(direc,1);
numRows = ceil(NFrames/numCols);
big_map_masks = double(zeros(H*numRows,W*numCols));
big_map_green = repmat(double(green),numRows,numCols);
r = 1;
c=1;
for n = 1:NFrames
    ims(:,:,n) = imread([direc(n).folder, '\',direc(n).name])-green;
    %[labs(:,:,n), regs(n)] = bwlabel(conv2(bwareaopen(conv2(ims(:,:,n),ones(5,5),'same'),900),ones(20,20), 'same'));
    [labs(:,:,n), regs(n)] = bwlabel(bwareaopen(ims(:,:,n),10));
    ax1 = axes;
    imagesc(green)
    colormap(ax1, 'gray');
    ax2 = axes;
    labsScaled = uint8(mat2gray(labs(:,:,n))*255);
    imagesc(ax2, labsScaled,'AlphaData',labs(:,:,n)>0);
    big_map_masks((r-1)*H+1:r*H,(c-1)*W+1:c*W) = labsScaled;
    colormap(ax2,jet);
    if(any(labsScaled>0))
        caxis(ax2,[min(nonzeros(labsScaled)) max(nonzeros(labsScaled))]);
    end
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    disp(['Frame: ', num2str(n), ' Groups: ',num2str(max(max(labs(:,:,n))))]);
    if(~exist([pathName,'\Colored\'],'dir'))
        mkdir([pathName,'\Colored']);
    end
    saveas(gcf,[pathName,'\Colored\',num2str(n)],'png');
    c = c+1;
    if(mod(n,numCols)==0)
        r = r+1;
        c = 1;
    end
end
figure();
plot(regs)
saveas(gcf,[pathName,'\Colored\Groups'],'fig');

small_map_green = imresize(big_map_green,0.25);
small_map_masks = imresize(big_map_masks,0.25);
figure();
ax_g = axes;
imagesc(small_map_green);
colormap(ax_g, 'gray');
ax_m = axes;
imagesc(ax_m, small_map_masks,'AlphaData',small_map_masks>0);
colormap(ax_m,jet);
caxis(ax_m,[min(nonzeros(small_map_masks)) max(nonzeros(small_map_masks))]);
ax_m.Visible = 'off';
linkprop([ax_g ax_m],'Position');
    
saveas(gcf, [pathName, '\Colored\All'], 'png');