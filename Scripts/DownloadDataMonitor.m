classdef DownloadDataMonitor < matlab.net.http.ProgressMonitor
    properties
        Value 
        Direction 
    end
    properties (Access = private)
        LastLength = 0;
        st;
    end
    methods
        function obj = DownloadDataMonitor()
            obj.Interval = 5;
            obj.st = tic;
            fprintf('\n');
        end
        function set.Value(obj, val)
            if obj.LastLength > 0
                fprintf(repmat('\b', 1, obj.LastLength));
            end
            msg = sprintf('Downloaded: %.2fGB', single(val)/1e9);
            if ~isempty(obj.Max)
                msg = [msg,sprintf('/%.2fGB %.1f%% (~%d s remaining)',single(obj.Max)/1e9,100*...
                    (single(val)/single(obj.Max)),fix((single(obj.Max)-single(val))/(single(val)/toc(obj.st))))];
            end
            fprintf('%s', msg);
            obj.LastLength = length(msg);
        end
    end
    
    methods (Access = public)
        function done(obj)
            fprintf("Done\n");
        end
    end
end
