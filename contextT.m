classdef contextT < context
    %CONTEXT Summary of this class goes here
    %   Detailed explanation goes here
    properties
        n1
        n2
    end
    
    methods
        function obj = contextT(p, q, n1, n2)
            %CONTEXT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 4
                n1 = [];
                n2 = [];
                if nargin < 2
                    p = [];
                    q = [];
                end
            end
            obj = obj@context(p, q);
            obj.n1 = n1;
            obj.n2 = n2;
            
        end
    end
end

