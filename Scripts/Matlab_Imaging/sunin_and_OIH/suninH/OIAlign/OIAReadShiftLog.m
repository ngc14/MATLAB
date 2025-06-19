function [filename, block, stim, frame, dx, dy, bestobj, reffilename, refstim, refframe, currentobj]=OIAReadShiftLog(filename)
%function shifts=OIAReadShiftLog(filename)

% Read a file contains alignment parameters
% Input: 'filename' a text file generaged by 'OIAWriteShiftLog.m'
% Output: 'shifts' is a structure contains following fields
%            'filename': 
%            'stim'
%            'frame': Above 3 parameters indicate which frame (be shifted)
%            'dx'
%            'dy': how much is shifted
%            'corr': correlation result or other objective function result
%            'reffilename'
%            'refstim'
%            'refframe': These 3 parameters indicate the base frame.            
%           example: shifts(3).filename(1

[filename, block, stim, frame, dx, dy, bestobj, reffilename, refstim, refframe, currentobj]= ...
    textread(filename, '%s %d %d %d %f %f %f %s %d %d %f');

return;


% any better way to do following cell assignment? c={filename, stim...} doesn't work
c=cell(size(stim,1), 9);

c(:,1)=filename;
c(:,2)=num2cell(stim);
c(:,3)=num2cell(frame);
c(:,4)=num2cell(dx);
c(:,5)=num2cell(dy);
c(:,6)=num2cell(corr);
c(:,7)=reffilename;
c(:,8)=num2cell(refstim);
c(:,9)=num2cell(refframe);
field={'filename', 'stim', 'frame', 'dx', 'dy', 'corr', 'reffilename', 'refstim', 'refframe'};
shifts=cell2struct(c, field,2);

return;

%doesn't work
%shifts=struct('filename', {filename}, 'stim', {stim}, 'frame', {frame}, 'dx', {dx}, 'dy', {dy}, 'corr', {corr}, 'reffilename', {reffilename}, 'refstim', {refstim}, 'refframe', {refframe});
