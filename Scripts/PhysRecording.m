classdef PhysRecording < BehavioralTask
    properties (GetAccess = 'public', SetAccess = 'private')
        binSize, sigmaSize, PSTHAlignments, secondsBeforePSTHAlignmentPoint,...
            secondsAfterPSTHAlignmentPoint;
    end
    properties (GetAccess = 'public', SetAccess = 'private', Dependent)
        sigma, bins;
    end
    properties (Access = 'private', Hidden)
        dPSTHBins = .01; % size of window to bin spikes (in seconds)
        dSigmaKernel = .15; % smoothing window in seconds
        % time (in seconds) before alignement point to create PSTH
        dSecondsBefore = -6;
        % time (in seconds) after alignement point to create PSTH
        dSecondsAfter = 5;
        % PSTH alignment points as segment names per condition
        dCondAlignments = containers.Map(["Extra Small Sphere", ...
            "Large Sphere", "Photocell", "Rest"],...
            {["StartReach", "StartHold", "StartWithdraw"], ...
            ["StartReach", "StartHold", "StartWithdraw"],...
            ["StartReach", "StartHold", "StartWithdraw"],...
            "GoSignal"});
    end
    
    methods (Static)
        function obj = PhysRecording(conditions,varargin)
            p = inputParser;
            addRequired(p,'conditions',@isstring);
            parse(p,conditions);
            obj@BehavioralTask(p.Results.conditions);
            validNumericParams = (@(x) isnumeric(x) & length(x)==1);
            addOptional(p,'binSize',obj.dPSTHBins,validNumericParams);
            addOptional(p,'sigma',obj.dSigmaKernel,validNumericParams);
            addOptional(p,'secondsBeforePSTH',obj.dSecondsBefore,validNumericParams);
            addOptional(p,'secondsAfterPSTH',obj.dSecondsAfter,validNumericParams);
            addOptional(p,'condAlignments',obj.dCondAlignments,(@(x) isobject(x) ...
                && isa(x,'containers.Map')));
            parse(p,p.Results.conditions,varargin{:});
            
            obj.binSize = p.Results.binSize;
            obj.sigmaSize = p.Results.sigma;
            obj.secondsBeforePSTHAlignmentPoint = p.Results.secondsBeforePSTH;
            obj.secondsAfterPSTHAlignmentPoint = p.Results.secondsAfterPSTH;
            obj.PSTHAlignments = containers.Map(p.Results.condAlignments.keys,...
                values(obj.dCondAlignments,p.Results.condAlignments.keys));
        end
    end
    methods
        function a = get.sigma(obj)
            a = fix(obj.sigmaSize./obj.binSize);
        end
        function b = get.bins(obj)
            b = obj.secondsBeforePSTHAlignmentPoint:obj.binSize:...
                (obj.secondsAfterPSTHAlignmentPoint-obj.binSize);
        end
    end
end