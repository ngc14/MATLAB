classdef MotorMapping
    properties (Constant)
        forelimbRep = ["Arm","Hand"];
        evokedResponseCategorization = containers.Map(["Digit", "Wrist", ...
            "Triceps", "Biceps","Elbow", "Forearm", "Shoulder","Face","Ear",...
            "Trunk", "Neck","No response"],...
            ["Hand", "Hand", "Arm", "Arm","Arm", "Arm", "Arm","Face", ...
            "Face", "Trunk", "Trunk", ""]);
        repColors = struct("Arm", [.85 .85 .85], "Hand", [.15 .15 .15],...
            "Forelimb",[.5 .5 .5], "Mixed", [225 150 165]./255,"Trunk",...
            [.5 .25 .1],"Face",[155 9 144]./255,"Axial",[102 0 36]./255);
%         repColors = struct("Reach", [1 0 0], "Grasp", [0 0 1], "Both", [1 0 1]);
        tileBuffer = 200;        
    end
    properties(Dependent)
        poolCircle
    end
    properties
        siteRadius;
    end
    methods
        function obj = MotorMapping(rd)
            if(~exist('rd', 'var'))
                rd = 42;
            end
            obj.siteRadius = rd;
        end
        function pc = get.poolCircle(obj)
            pc = fspecial('disk',obj.siteRadius);
        end
    end
end
