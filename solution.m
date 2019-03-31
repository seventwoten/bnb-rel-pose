classdef solution
    %SOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        centre   % block centre
        sigma    % block half-length 
        score    % block upper bound
        edges    % edges matrix
        
    end
    
    methods
        function obj = solution(c,s,m,e)
            %SOLUTION Construct an instance of this class
            %   Detailed explanation goes here
            obj.centre = c;
            obj.sigma  = s;
            obj.score  = m;
            obj.edges  = e;
            
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

