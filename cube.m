classdef cube < block
    %CUBE 3D block
    methods
        function obj = cube(c, s, p)
            %CUBE Construct an instance of this class
            if nargin == 0
                c = [];
                s = [];
                p = [];
            end
            obj = obj@block(c, s, p);
        end
    
        function [subblocks] = subdivide(obj)
            %SUBDIVIDE Subdivide block and return 8 subblocks
            shifts = [-1, -1, -1 ;
                      -1, -1,  1 ;
                      -1,  1, -1 ;
                      -1,  1,  1 ;
                       1, -1, -1 ;
                       1, -1,  1 ;
                       1,  1, -1 ;
                       1,  1,  1  ];
            
            sigma_new  = obj.sigma/2;
            centre_new = obj.centre + sigma_new * shifts;
            
            for c = size(centre_new, 1):-1:1
                subblocks(c) = cube(centre_new(c,:), sigma_new, obj.UB);
            end
            
        end
    end
end

