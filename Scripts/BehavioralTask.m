classdef BehavioralTask
    properties (Hidden, Constant)
        dCondSegs = containers.Map(...
            ["Extra Small Sphere", "Large Sphere","Photocell","Rest"], ...
            {["StartTrial","GoSignal","StartReach","StartGrasp","StartLift","StartHold",...
            "StartWithdraw","StartReplaceHold","StartReplaceSuccess","StartReward","EndTrial"],...
            ["StartTrial","GoSignal","StartReach","StartGrasp","StartLift","StartHold",...
            "StartWithdraw","StartReplaceHold","StartReplaceSuccess","StartReward","EndTrial"],...
            ["StartTrial","GoSignal","StartReach","StartGrasp","StartHold","StartWithdraw",...
            "StartReplaceHold","StartReplaceSuccess","StartReward","EndTrial"],...
            ["StartTrial","GoSignal","StartReplaceHold","StartReplaceSuccess","StartReward","EndTrial"]});
    end
    properties (GetAccess = 'public', SetAccess = 'private')
        condNames,segNames,condAbbrev,condSegMap;
    end
    methods
        function obj = BehavioralTask(varargin)
            if nargin == 2 && isstring(varargin{1}) && isstring(varargin{2})...
                    && length(varargin{1})==length(varargin{2})
                obj.condNames = varargin{1};
                obj.segNames = varargin{2};
            else
                if nargin == 0
                    obj.condNames = obj.dCondSegs.keys();
                elseif nargin == 1 && isstring(varargin{1})
                    obj.condNames = varargin{1};
                end
                obj.segNames = values(obj.dCondSegs,cellstr(obj.condNames));
            end
            obj.condSegMap = containers.Map(obj.condNames,obj.segNames);
            obj.condAbbrev = containers.Map(obj.condNames,cellfun(@(cc) ...
                cc([1,(regexp(cc,'\s'))+1]),obj.condNames, 'UniformOutput', false));
        end
    end
end