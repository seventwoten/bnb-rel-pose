classdef block
    %BLOCK Represents a 2D/3D region for use in branch-and-bound
    
    properties
        centre   % block centre
        sigma    % block half-length 
        parentUB % parent block upper bound
        
        UB       % block upper bound
        LB       % block lower bound
    end
    
    methods
        function obj = block(c, s, p)
            %BLOCK Construct an instance of this class
            if nargin ~= 0
                obj.centre = c;
                obj.sigma = s;
                obj.parentUB = p;
            end
        end
        
    end
end

