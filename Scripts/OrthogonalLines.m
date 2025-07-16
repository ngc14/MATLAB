classdef OrthogonalLines
    properties
        ArmTrunkCS;
        ArmFaceCS;
    end
    properties (Dependent)
        RCLine
        MLLine
    end
    methods
        function mn = OrthogonalLines(monkey)
            if(strcmp(monkey,"Gilligan"))
                mn.ArmTrunkCS=[30,70];
                mn.ArmFaceCS=[75,536];
            elseif(strcmp(monkey,"Skipper"))
                mn.ArmTrunkCS = [266,226];
                mn.ArmFaceCS = [595,637];
            end
        end
        function rcVals = get.RCLine(obj)
            rcVals = [obj.ArmTrunkCS(1),obj.ArmTrunkCS(end),obj.ArmFaceCS(1),obj.ArmFaceCS(end)];
        end
        function mlVals = get.MLLine(obj)
            findFit1 = polyfit([obj.ArmTrunkCS(1),obj.ArmFaceCS(1)],[obj.ArmTrunkCS(end), obj.ArmFaceCS(end)],1);
            newSlope = -1/findFit1(1);
            mlVals = [obj.ArmTrunkCS(1),obj.ArmTrunkCS(end),obj.ArmFaceCS(1),...
                newSlope*(obj.ArmFaceCS(end)-obj.ArmTrunkCS(1)) + obj.ArmTrunkCS(2)];
        end
    end
end