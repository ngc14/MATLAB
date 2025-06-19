function IMManualClip()

clear all;
% manually adjust clip range
% "up/down" keys to adjust mean
% "left/right" keys to adjust clip width (stand deviation), right key to increase width
% "space" key to save
% "esc" key to quite (no save)

filename='g0.ivf';
step=0.1;       % each key push will adjust std*step 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map=OIReadIVF(filename);
mapmean=mean(map(:));
mapstd=std(map(:));
step=mapstd*step;
clipmin=mapmean-mapstd;
clipmax=mapmean+mapstd;
if ~isempty(dir('clipvalues.txt'))
    temp=textread('clipvalues.txt', '%f');
    clipmin=temp(1);
    clipmax=temp(2);
end

BKey = KbName('b');	% 
rightKey = KbName('right');	% 
leftKey = KbName('left');	% 
upKey = KbName('up');
downKey = KbName('down');
shiftKey = kbName('left_shift');   % left shift
ctrlKey = kbName('left_control'); % left ctrl
altKey = kbName('alt');   % left alt
pageupKey = kbName('pageup');
pagedownKey = kbName('pagedown');
spaceKey = kbName('space');
escapeKey = KbName('esc');	% quit key

fprintf('<-, -> or Esc\n');
n=0;
while 1
    image(norm_to_uint8b(map, clipmin, clipmax))
    axis equal;
    colormap(gray(256));
    drawnow;
    % pause(0.01);
    [keyIsDown,secs,keyCode] = KbCheck;
    if keyIsDown
        if keyCode(rightKey) 
            clipmin=clipmin-step;
            clipmax=clipmax+step;
        elseif keyCode(leftKey)
            clipmin=clipmin+step;
            clipmax=clipmax-step;
        elseif keyCode(upKey)
            clipmin=clipmin+step;
            clipmax=clipmax+step;
        elseif keyCode(downKey)
            clipmin=clipmin-step;
            clipmax=clipmax-step;
        elseif keyCode(spaceKey)
            n=n+1;
            imwrite(norm_to_uint8b(map, clipmin, clipmax), ['final',num2str(n), '.bmp']);
            fid=fopen('clipvalues.txt', 'w');
            fprintf(fid,  '%f\t%f', clipmin, clipmax);
            fclose(fid);
            fprintf('current map saved as "%d.bmp" and setting in "clipvalues.txt"\r', n);
            
        end
        if keyCode(escapeKey)      % 'Esc' key for Exit
            break;
        end
    end % if KeyIsDown loop
end % while loop
return;      
