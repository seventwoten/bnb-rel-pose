classdef sphericalPatch < block
    %SPHERICALPATCH 2D patch on spherical surface 
    properties
        angleMat1 % angles between (t, n1)
        angleMat2 % angles between (t, n2)
        centre_xyz % cartesian representation of patch centre
    end
    
    methods
        function obj = sphericalPatch(c, s, p)
            %SPHERICALPATCH Construct an instance of this class
            if nargin == 0
                c = [];
                s = [];
                p = [];
            end
            obj = obj@block(c, s, p);
            
            if nargin == 3
                obj.thres = acos(cos(s)^2);
                obj.centre_xyz = obj.spherical2Cartesian(1, c(1), c(2));
            end
        end
        
        function [subblocks] = subdivide(obj)
            %SUBDIVIDE Subdivide block and return 8 subblocks
            shifts = [-1, -1 ;
                      -1,  1 ;
                       1, -1 ;
                       1,  1  ];
            
            sigma_new  = obj.sigma/2;
            centre_new = obj.centre + sigma_new * shifts;
            
            for c = size(centre_new, 1):-1:1
                subblocks(c) = sphericalPatch(centre_new(c,:), sigma_new, obj.UB);
            end
            
        end
    end
    
    methods (Static)
        function [xyz] = spherical2Cartesian(r, theta, phi)
            sp = sin(phi);
            cp = cos(phi);
            st = sin(theta);
            ct = cos(theta);

            x = r * sp * ct;
            y = r * sp * st;
            z = r * cp;
            xyz = [x, y, z];
        end
    end
end


