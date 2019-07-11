classdef tPatch < block
    %TPATCH 2D patch on spherical surface 
    properties
        centre_xyz % cartesian representation of patch centre
    end
    
    methods
        function obj = tPatch(c, s)
            %TPATCH Construct an instance of this class
            if nargin == 0
                c = [];
                s = [];
            end
            obj = obj@block(c, s);
            
            % Compute thres and centre_xyz if inputs are not []
            if nargin == 2
                obj.thres = acos(cos(s)^2);
                obj.centre_xyz = obj.spherical2Cartesian(1, c(:,1), c(:,2));
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
                subblocks(c) = tPatch(centre_new(c,:), sigma_new);
            end
            
        end
    end
    
    methods (Static)
        function [xyz] = spherical2Cartesian(r, theta, phi)
            %SPHERICAL2CARTESIAN Convert spherical coordinates to Cartesian
            %   Can convert multiple vectors by passing in multiple r, 
            %   theta and phi as column vectors
            sp = sin(phi);
            cp = cos(phi);
            st = sin(theta);
            ct = cos(theta);

            x = r .* sp .* ct;
            y = r .* sp .* st;
            z = r .* cp;
            xyz = [x, y, z];
        end
    end
end


