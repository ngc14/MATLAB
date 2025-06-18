%ASSUMPTIONS
% -MOTOR MAPPING INFO
%   - \dirPath\Mapping\*_MMSites.xlsx : file that stores all recording and 
%   and motor mapped sites for monkey. Header of file contains:
%   Recording (Y/N), Thresholds, Evoked Movement, X site location, Y site location
%
% -M1+PREMOTOR BORDER   
%   - \dirPath\Mapping\*_MMRGB.png : image of FOV with M1/premotor border 
%   drawn in yellow
%
% -TTESTS  
%   - \dirPath\Mapping\*tTest.png : images of pixels (red) that were significant
%   from imaging analysis 
%
% -RECORDING INFO
%   - \dirPath\*_Session.xlsm : file that stores session information. Header
%   of file contains: Date, Site #, and Domain (brain area) of the recording 
%
%RETURNS
% -infoTable
%   - table of recording sites in M1 (validated by M1+PREMOTOR BORDER) with 
%   site #, location, date, evoked movement response at neartest MM site, 
%   minimum current threshold at nearest MM site
%
% -masks
%   - FOV blood vessel mask, M1, PMd and PMv mask
%
% -activityMaps (optional)
%   - images in the same FOV that show condition specific activity modulation
%
function [infoTable,masks,activityMaps] = getMonkeyInfo(dirFolder, monkey,domains,recordingsOnly)
dirPath = strcat(dirFolder,monkey,"\");
fileName = monkey;
%%% vessel mask
vMask =  im2double(im2gray(imread(strcat(dirPath, "Mapping\clean_reference_mask.png"),"BackgroundColor","none")));
%%% all mapped sites and x,y locations
rawRef=readcell(strcat(dirPath,"Mapping\",fileName,"_MM_Sites.xlsx"));
% file header info
heading = rawRef(1,:);
rawRef(1,:) = [];
[~,recordingInd] = find(strcmp(heading, "Recording"));
[~,siteInd] = find(strcmp(heading, "Site #"));
[~,threshInd] = find(strcmp(heading, "Threshold(s)"));
[~,erInd] = find(strncmp(heading, "Evoked",length("Evoked")));
[~,xInd] = find(strcmp(heading, "x"));
[~,yInd] = find(strcmp(heading, "y"));
% get recording sites, not MM sites
if(recordingsOnly)
    validSiteInds = cellfun(@(a,b,c,d) strcmp(a,"Yes") & isnumeric(b) & ~ismissing(c)...
        & ~ismissing(d),rawRef(:,recordingInd),rawRef(:,siteInd),rawRef(:,xInd),rawRef(:,yInd));
else
    validSiteInds = cellfun(@(a,b,c) isnumeric(a) & ~ismissing(b) & ~ismissing(c), ...
        rawRef(:,siteInd),rawRef(:,xInd),rawRef(:,yInd));
end
% if(strcmp(monkey,"Skipper"))
%     validSiteInds = validSiteInds & [rawRef{:,siteInd}]'<247;
% elseif(strcmp(monkey,"Gilligan"))
%     validSiteInds = validSiteInds & [rawRef{:,siteInd}]'<290 & [rawRef{:,siteInd}]'>6;
% end
% store locaiton and site number for each site
infoTable = array2table(cellfun(@str2double, string(rawRef(validSiteInds,...
    [siteInd,xInd,yInd]))),'VariableNames', ["Site", "x", "y"]);
% mapped sites to determine identity of recorded sites
mappedSiteInds = find(~cellfun(@(e) (length(e)==1 && isnan(e)) || ...
    strcmp(e, "No response"), rawRef(:,erInd)));
allLocs = [cellfun(@double, rawRef(mappedSiteInds,xInd), 'UniformOutput', true),...
    cellfun(@double, rawRef(mappedSiteInds,yInd), 'UniformOutput', true)];
validSiteInds = find(validSiteInds);
[allReps, allThreshs] = deal(cell(length(validSiteInds),1));
for r = 1:length(validSiteInds)
    currLoc = [double(rawRef{validSiteInds(r),xInd}), ...
        double(rawRef{validSiteInds(r),yInd})];
    % find mapped site closest to recorded site
    [~,minMappedInd] = min(sqrt(sum((allLocs-currLoc).^2,2)));
    % map evoked response names and current thresholds of closest mapped
    % site to current recorded site
    currERs = regexp(rawRef{mappedSiteInds(minMappedInd),erInd},',','split');
    thresh = rawRef{mappedSiteInds(minMappedInd),threshInd};
    erLabels = MotorMapping.evokedResponseCategorization.keys();
    mappedERs = cellfun(@(e) find(cellfun(@(ab) contains(strip(e),ab), ...
        erLabels),1),currERs);
    allReps{r} = string(values(MotorMapping.evokedResponseCategorization,...
        erLabels(mappedERs)));
    if(length(thresh)<length(mappedERs))
        disp(rawRef{mappedSiteInds(minMappedInd),1});
    end
    if(isnumeric(thresh))
        allThreshs{r} = thresh;
    else
        allThreshs{r} = cellfun(@(t) str2double(strip(t)),...
            regexp(thresh, ',', 'split'));
    end
end
infoTable = addvars(infoTable,allReps, allThreshs,repmat(monkey,...
    length(allThreshs),1),'NewVariableNames',["SiteRep", "Thresh","Monkey"]);
clear rawNum rawRef
%%%%%%%%%% date and session information (double check M1 mask)%%%%%%%%%%
opts = detectImportOptions(strcat(dirPath,fileName,"Session.xlsm"));
opts.VariableNamingRule = 'preserve';
opts.Sheet=1;
opts.DataRange='B1';  
opts.VariableNamesRange='' ; 
opts=setvaropts(opts,'EmptyFieldRule','auto');
dateToPNValue = readcell(strcat(dirPath,fileName,"Session.xlsm"),opts);
[~, dateInd] = find(strcmp(dateToPNValue, "Date"));
[~,domainInd] = find(strcmp(dateToPNValue, "Domain"));
[~,PNInd] = find(strcmp(dateToPNValue, "Site #"));
[~,notesInd] = find(strcmp(dateToPNValue, "NOTES:"));

% clean up, remove: duplicate test sites (close/far), non-motor recordings,
% improper logging (no site #)
validInds = cellfun(@(a,b,c,d) all(~ismissing(a)) && all(~ismissing(b)) &&  ...
    any(contains(b, domains)) && all(~ismissing(c)) && ~contains(string(d),"Single"),...
    dateToPNValue(:,dateInd),dateToPNValue(:,domainInd),dateToPNValue(:,PNInd),dateToPNValue(:,notesInd)); 

% match session numbers from _MM_Sites.xlsx sheet to session.xlsx sheet
[~,viPN,vI] = intersect(cellfun(@str2double, string(dateToPNValue(...
    validInds,PNInd))),infoTable.Site);
% infoTable = infoTable(viPN,:);
dateToPNValue = dateToPNValue(validInds,:);
infoTable = [infoTable table(cell(height(infoTable),1),...
    cell(height(infoTable),1),zeros(height(infoTable),1),'VariableNames',...
    {'Date', 'Domain', 'Site2'})];
% assign session number's date, add domain, duplicate session number from
% session.xlsx sheet
infoTable(vI,'Date') = table(dateToPNValue(viPN,dateInd));
infoTable(vI,'Domain') = table(dateToPNValue(viPN,domainInd));
infoTable(vI,'Site2')= dateToPNValue(viPN,PNInd);
clear dateToPNValue
%%%%%%% confirm M1 sites based on assigned yellow border border
MM = double(imread(strcat(dirPath,"Mapping\MMRGB.png")));
M1border = bwareafilt(MM(:,:,1)>125 & MM(:,:,2)>125 & MM(:,:,3)<125, 1);
PMborder = bwareafilt(MM(:,:,1)>100 & MM(:,:,2)<125 & MM(:,:,3)>100, 1);
% find border
[row,col] = find(M1border);
pointsX = unique(row);
pointsY = arrayfun(@(a) max(col(row==a)),pointsX);
pointsX(end+1:end+3) = [size(MM,1),size(MM,1),0];
pointsY(end+1) = pointsY(end);
pointsY(end+1:end+2) = [0,0];
% create M1/PM mask
M1mask = double(poly2mask(pointsY,pointsX, size(MM,1),size(MM,2)));
[rowPM,colPM] = find(PMborder);
pointsXPM = unique(rowPM);
pointsYPM = arrayfun(@(a) max(colPM(rowPM==a)),pointsXPM);
pointsXPM(end+1:end+2) = [size(MM,1) size(MM,1)];
pointsYPM(end+1:end+2) = [size(MM,2) 0];
PMVmask = double(~M1mask & poly2mask(pointsYPM,pointsXPM,size(MM,1),size(MM,2)));
PMDmask = ~M1mask & ~PMVmask;
masks = {vMask,M1mask,PMDmask,PMVmask};
% validate M1 sites based on location
infoTable{find(arrayfun(@(a,b) M1mask(b,a), infoTable.x,infoTable.y)),'Domain'} = {"M1"};
infoTable{find(arrayfun(@(a,b) PMDmask(b,a), infoTable.x,infoTable.y)),'Domain'} = {"PMd"};
infoTable{find(arrayfun(@(a,b) PMVmask(b,a), infoTable.x,infoTable.y)),'Domain'} = {"PMv"};
infoTable(find(cellfun(@(s) ~any(strcmp(s,domains)), infoTable.Domain)),:) = [];
%infoTable([infoTable.Site]~=[infoTable.Site2],:) = [];
%%% activity maps from imaging
allActivityFiles = dir(strcat(dirPath,"Mapping\*tTest*.bmp"));
allMaps = arrayfun(@(a) imresize(imread(strcat(a.folder, "\", a.name)),...
    [size(MM,1),size(MM,2)]), allActivityFiles, 'UniformOutput', false);
% activity is red
activityMaps = containers.Map();
for ma = 1:length(allMaps)
    activityMaps(strcat(monkey,extractBefore(allActivityFiles(ma).name,"_tTest")))=...
        allMaps{ma}>0;
end
end