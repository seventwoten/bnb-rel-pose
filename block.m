classdef block
    %BLOCK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        centre   % block centre
        sigma    % block half-length 
        parentUB % parent block upper bound
        
    end
    
    methods
        function obj = block(c, s, p)
            %BLOCK Construct an instance of this class
            %   Detailed explanation goes here
            obj.centre = c; 
            obj.sigma = s;  
            obj.parentUB = p;
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

