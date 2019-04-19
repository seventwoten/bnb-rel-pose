classdef block
    %BLOCK Represents a 2D/3D region for use in branch-and-bound
    
    properties
        centre   % block centre
        sigma    % block half-length 
        parentUB % parent block upper bound
        
        dim      % dimension (3 = cube, 2 = sphere surface patch)
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
                obj.dim = numel(c);
                assert(obj.dim == 2 | obj.dim == 3, 'Block centre is not 2D or 3D.')
            end
        end
        
        function [subblocks] = subdivide(obj)
            %SUBDIVIDE Subdivide block and return subblocks
            
            if obj.dim == 3
                shifts = [-1, -1, -1 ;
                          -1, -1,  1 ;
                          -1,  1, -1 ;
                          -1,  1,  1 ;
                           1, -1, -1 ;
                           1, -1,  1 ;
                           1,  1, -1 ;
                           1,  1,  1  ];

            else % obj.dim == 2
                shifts = [-1, -1 ;
                          -1,  1 ;
                           1, -1 ;
                           1,  1  ];
            end
            
            sigma_new  = obj.sigma/2;
            centre_new = obj.centre + sigma_new * shifts;
            
            for c = size(centre_new, 1):-1:1
                subblocks(c) = block(centre_new(c,:), sigma_new, obj.UB);
            end
            
        end
    end
end

