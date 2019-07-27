classdef tPatch < block
    %TPATCH 2D patch on spherical surface 
    properties
        centre_xyz % cartesian representation of patch centre
        R_tilt     % Rotation around origin to reposition the patch on the 
                   % sphere surface, used if patch is not aligned to long/lat lines
    end
    
    methods
        function obj = tPatch(c, s, R)
            %TPATCH Construct an instance of this class
            if nargin == 0
                c = [];
                s = [];
            end
            obj = obj@block(c, s);
            
            % Compute thres and centre_xyz if inputs are not []
            if ~isempty(s)
                assert(s <= pi/2, 'Sigma (patch half-length) must not exceed pi/2');
                obj.thres = acos(cos(s)^2);
            end
            if ~isempty(c)
                obj.centre_xyz = obj.spherical2Cartesian(1, c(:,1), c(:,2));
            end
            
            % If R_tilt is given: Rotate patch, update centre and centre_xyz
            if exist('R','var') && ~isempty(R)
                obj.R_tilt = R;
                obj = obj.rotate(R);
            end
        end
        
        function [subblocks] = subdivide(obj)
            %SUBDIVIDE Subdivide block and return 8 subblocks
            if ~isempty(obj.R_tilt)
                obj = obj.rotate(obj.R_tilt');
            end
            shifts = [-1, -1 ;
                      -1,  1 ;
                       1, -1 ;
                       1,  1  ];
            
            sigma_new  = obj.sigma/2;
            centre_new = obj.centre + sigma_new * shifts;
            
            for c = size(centre_new, 1):-1:1
                subblocks(c) = tPatch(centre_new(c,:), sigma_new, obj.R_tilt);
            end
            
        end
        
        function [obj] = rotate(obj, R)
            %SUBDIVIDE Rotate the patch by R, and update centre_xyz, centre
            obj.centre_xyz = obj.centre_xyz * R';
            rtp = obj.cartesian2Spherical(obj.centre_xyz(:,1), obj.centre_xyz(:,2), obj.centre_xyz(:,3));
            obj.centre = [rtp(:,2), rtp(:,3)];
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
        
        function [rtp] = cartesian2Spherical(x, y, z)
            %CARTESIAN2SPHERICAL Convert Cartesian to spherical coordinates
            %   Can convert multiple vectors by passing in multiple x, y, z
            %   as column vectors
            
            r     = sqrt(x.^2 + y.^2 + z.^2);
            theta = atan2(y,x);
            phi   = atan2(sqrt(x.^2 + y.^2),z);
            rtp = [r, theta, phi];
        end
    end
end


