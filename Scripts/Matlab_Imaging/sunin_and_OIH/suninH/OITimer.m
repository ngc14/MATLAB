function [ctime, timertrigger]=OITimer(tnow, t0, timertrigger, Nblock, NFramePerBlock)

%% Just for estiming time to do processing
%% It can be used in a loop without any change of it's format, i.e.:
%% 		[ctime, timertrigger]=OITimer(clock, ctime, timertrigger, Nblock, NFramePerBlock)
%% first time need set timertrigger=1 
%% 2nd time make no change, when enter this sub-routine, it print out the time it used 
if timertrigger==1
    ctime=tnow;
    timertrigger=2;
    return;
elseif timertrigger==2
    frametime=etime(tnow, t0);
    ctime=tnow; %no use
    timertrigger=0;
    fprintf('\rTime for one frame: %7.4f secs', frametime);
    fprintf('\rTime for one block: %5.2f mins', frametime*NFramePerBlock/60);
    fprintf('\rTime for all (%d) blocks: %5.2f hours\r', Nblock, frametime*NFramePerBlock*Nblock/3600);
    return
end