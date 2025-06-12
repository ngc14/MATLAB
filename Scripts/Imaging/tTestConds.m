clear all; close all;
monkey = "Gilligan";
rootFolder = "E:";
conds = ["[ExtraSmallSphere]","[LargeSphere]", "[Photocell]"];
condComb = nchoosek(1:length(conds),2);
HPk = 250;
LPk = 5;
pval = 0.001;
cleanArea = 10;
saveFold = "C:\Lab\"+monkey+"\tTestConds\";
if(strcmp(monkey,"Gilligan"))
    dates = ["12_05_2018", "12_07_2018","12_09_2018", "12_10_2018", "12_12_2018", "12_13_2018","12_14_2018","01_07_2019"];
    runs = [0, 0, 0, 1, 0, 0, 1, 0];
    frames = 49:53;
elseif(strcmp(monkey,"Skipper"))
    dates = ["10_30_2020","11_09_2020","11_23_2020", "11_24_2020", "11_25_2020", "11_27_2020", "11_30_2020","12_01_2020","12_02_2020"];
    runs = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    frames = 48:52;
end
frames = 32:70;
refMask = imbinarize(imread(rootFolder+"\"+monkey+"\Mapping\clean_mask_filled",'bmp'));
refMask = refMask(:,:,1);
%%
ims = cell(1,length(conds));
for c = 1:length(conds)
    currTrials= loadAllTrials(rootFolder,monkey,dates,runs,LPk,HPk,conds(c),frames);
    ims{c} = tall(cat(3,currTrials{:}));
    disp(toc);
end
%%
for c = 1:length(condComb)
    cmpIm = ones(size(refMask,[1 2])).*5;
    if(~exist(saveFold,'dir'))
        mkdir(saveFold);
    end
    fileName = saveFold + conds(condComb(c,1))+ "_"+ conds(condComb(c,2)) +"_"+ ...
        num2str(frames(1)) + "_" + num2str(frames(end)) + "_" + num2str(pval);
    im1 = gather(ims{condComb(c,2)});
    im2 = gather(ims{condComb(c,1)});
    [~,imP1] = ttest2(im1,im2,'alpha',0.05,'tail','right','dim',3);
    [~,imP2] = ttest2(im2,im1,'alpha',0.05,'tail','right','dim',3);
    cmpIm(bwareaopen(imfill(imP1 <pval,'holes'),cleanArea)) = 2;
    cmpIm(bwareaopen(imfill(imP2 <pval,'holes'),cleanArea)) = 3;
    cmpIm(bwareaopen(imfill(imP1<pval & imP2<pval,'holes'),cleanArea)) = 4;
    cmpIm(~refMask) = 1;
    figure();
    imagesc(cmpIm); colormap([.5 .5 .5; 1 0 0; 0 1 0; 0 0 1; 1 1 1]);
    clim([1 5])
    title([conds{condComb(c,1)}, '(R) > ', conds{condComb(c,2)}, '(G)'])
    g = gca;
    imwrite(g.Children.CData,colormap,fileName + '.png');
    exportgraphics(g,fileName + '.eps','ContentType','vector');
end