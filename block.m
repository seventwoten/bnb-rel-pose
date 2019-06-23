classdef block
    %BLOCK Represents a 2D/3D region for use in branch-and-bound
    
    properties
        centre   % block centre
        sigma    % block half-length 
        
        UB       % block upper bound
        LB       % block lower bound
        
        thres       % threshold for comparison with stopping threshold
        edges_stop  % edges matrix of possible matches at stopping threshold
    end
    
    methods
        function obj = block(c, s)
            %BLOCK Construct an instance of this class
            if nargin ~= 0
                obj.centre = c;
                obj.sigma = s;
            end
        end
    end
    
end

