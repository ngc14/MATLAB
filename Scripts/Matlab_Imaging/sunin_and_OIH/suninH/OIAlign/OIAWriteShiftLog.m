function OIAWriteShiftLog(fid, filename, block, stim, frame, dx, dy, bestobj, reffilename, refstim, refframe, currentobj)

% OIAWriteShiftLog(fid, filename, block, stim, frame, dx, dy, corr, reffilename, refstim, refframe)
% Write an entry in shift log.

fprintf(fid, '%s \t%d \t%d \t%d \t%3.1f \t%3.1f \t%6.6f \t%s \t%d \t%d  \t%f \r\n', filename, block, stim, frame, dx, dy, bestobj, reffilename, refstim, refframe, currentobj);

return
