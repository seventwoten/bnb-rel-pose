classdef context
    %CONTEXT Summary of this class goes here
    %   Detailed explanation goes here
    properties
        p
        q
    end
    
    methods
        function obj = context(p, q)
            %CONTEXT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                p = [];
                q = [];
            end
            
            obj.p = p;
            obj.q = q;
            
        end
    end
end

