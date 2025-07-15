close all;
clear all;
%%WINDOWS MAGNIFICATION MUST BE SET TO 100%
oldVidPath = 'S:\Lab\ngc14\Track_Videos\ISOI_Behave_Vid_Sessions.avi';

behaveFR = 120;
ranges = {[4.5, 5.0], [5.0, 5.4], [5.4, 5.8], [5.8, 7.1]};
domainsPerRange = {{'M1_CS', 'M1/PMd'}, {'M1_CS', 'M1/PMd', 'PMd_L'}, {'M1_CS', 'M1/PMd', 'PMd_L', 'PMd_M'},...
    {'M1/PMd', 'PMd_L', 'PMd_M', 'M1_CS_R'}};
domainPositions =  containers.Map({'M1_CS', 'M1/PMd', 'PMd_L', 'PMd_M', 'M1_CS_R'}, ...
{{[.675, .674], [.485, .5]}, {[.86, .893], [.67, .57]}, {[.71,.73], [.73, .63]}, {[.87, .90], [.72, .728]},...
{[.674, .676], [.485, .51]}});
screensize = get( groot, 'Screensize' );
oldVid = VideoReader(oldVidPath);
newVid = VideoWriter([oldVidPath(1:end-4), '_Annotated.avi']);
newVid.FrameRate = 60;
open(newVid);
count = 1;

while(oldVid.hasFrame)
    currFrame = im2double(oldVid.readFrame);
    vidTime = round(count/behaveFR,1);
    validRange = cellfun(@(a) vidTime >= a(1) && vidTime < a(2), ranges);
    
    if(find(validRange))
        f1 = figure('units','pixel','outerposition',[100, 100, size(currFrame,2), size(currFrame,1)]);
        set(f1,'Renderer','ZBuffer');
        iptsetpref('ImshowBorder','tight');
        imshow(currFrame, 'Border', 'tight');
               
        annPositions = cellfun(@(a) domainPositions(a), domainsPerRange{validRange}, 'UniformOutput', false);
        
        for n = 1:length(annPositions)
            currAnnotation = annPositions{n};
            annotation(gcf, 'arrow', annPositions{n}{1}, annPositions{n}{2},'Color', 'black', ...
                'HeadStyle', 'cback2', 'LineStyle', 'none', 'Units', 'pixels', ...
                'HeadLength', 30, 'HeadWidth', 25);
        end
        
        pause(0.2);
        grabFrame =getframe(gcf);
        currFrame = im2double(grabFrame.cdata);
        close all;
    end
    
    writeVideo(newVid,im2uint8(currFrame));
    count = count + 1;
end
close(newVid);