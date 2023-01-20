classdef MagneticField
    properties
        direction = [0,0];
        strength = 0;
        minAngle = 30;
        maxAngle = 150;
    end
    
    methods
        function obj = MagneticField(direction,strength,min,max)
            obj.direction = direction;
            obj.strength = strength;
            obj.minAngle = min;
            obj.maxAngle = max;
        end
    end
end