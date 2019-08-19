classdef wedge < handle 
    %WEDGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n1      % 1st normal (unit row vector)
        n2      % 2nd normal (unit row vector)
        centre  
        corner  
        width
    end
    
    methods
        function obj = wedge(normal1,normal2)
            %WEDGE Construct an instance of this class
            assert(obj.isValid(normal1, normal2), "Normals must be non-zero 3D vectors.")
            obj.n1 = normal1(:)';
            obj.n1 = obj.n1 ./ sqrt(sum(obj.n1.^2));
            obj.n2 = normal2(:)';
            obj.n2 = obj.n2 ./ sqrt(sum(obj.n2.^2));
        end
        
        function setCentre(obj)
            obj.centre = (obj.n1 + obj.n2) / 2; 
            obj.centre = obj.centre ./ sqrt(sum(obj.centre.^2));
        end
        
        function setCorner(obj)
            %SETCORNER Sets one wedge corner (the other is simply -corner)
            if isempty(obj.centre)
                obj.setCentre();
            end
            obj.corner = cross(obj.centre, obj.n1); 
            obj.corner = obj.corner ./ sqrt(sum(obj.corner.^2));
        end
        
        function setWidth(obj)
            obj.width  = abs(acos(obj.n1 * -obj.n2'));
        end
        
        function centre = getCentre(obj)
            if isempty(obj.centre)
                obj.setCentre();
            end
            centre = obj.centre;
        end
        
        function corner = getCorner(obj)
            if isempty(obj.corner)
                setCorner(obj);
            end
            corner = obj.corner;
        end
        
        function width = getWidth(obj)
            if isempty(obj.width)
                setWidth(obj);
            end
            width = obj.width;
        end
        
        function t_list = wedge2tPatches(obj, min_size, join_patches)
            %SETCONTEXT Convert a single wedge into a list of tPatches.
            %   Align wedge to z-axis, divide into patches, and rotate back.
            
            if ~exist('join_patches', 'var') || isempty(join_patches)
                join_patches = true;
            end
            % Find wedge_centre coords after aligning wedge corner to z-axis
            R_align = RCube.align2R(obj.getCorner(), [0,0,1]);
            wedge_centre_aln = R_align *  obj.getCentre()';
            
            % Find number and size of patches along wedge width
            min_size = max(pi/32, min_size);
            num_theta      = ceil( obj.getWidth() / (2*min_size) );
            patch_half_len = min_size;
            
            if join_patches
                while mod(num_theta, 2) == 0
                    num_theta      = num_theta  / 2;
                    patch_half_len = patch_half_len * 2;
                end
            end
            % Find number of patches along wedge length (even number)
            num_phi = round( pi / (2*patch_half_len) );
            
            % Initialise tPatches aligned to z-axis, then rotate back
            rtp = tPatch.cartesian2Spherical(wedge_centre_aln(1),wedge_centre_aln(2),wedge_centre_aln(3)); 
            thetas = rtp(2) + (-floor(num_theta/2):floor(num_theta/2)) * patch_half_len * 2;
            phis   = rtp(3) + ((1 : num_phi) - floor(num_phi/2) - 1) * patch_half_len * 2 + patch_half_len;
            
            [t,p] = meshgrid(thetas,phis);
            t_centres = [t(:), p(:)];
            t_list = tPatchList(t_centres, patch_half_len, R_align');
        end
        
    end
    methods (Static)
        function result = isValid(normal1, normal2)
            %ISVALID Checks if wedge normals are zero vectors
            is_3d = numel(normal1) == 3 && numel(normal2) == 3;
            not_zero = ~all(normal1 == 0) && ~all(normal2 == 0);
            result = is_3d && not_zero;
        end
        
        function int_points = wedges2int(w1, w2)
            % WEDGES2INTERSECTION Return intersection vertices and its mean
            
            % Find 8 intersections between 2 sets of 2 great circles
            intn = zeros(8, 3);
            intn(1,:) = cross(w1.n1, w2.n1);
            intn(2,:) = cross(w1.n1, w2.n2);
            intn(3,:) = cross(w1.n2, w2.n1);
            intn(4,:) = cross(w1.n2, w2.n2);
            intn(5:8, :) = -intn(1:4, :);
            intn = intn ./ sqrt(sum(intn.^2, 2));            

            % Find wedge centres and corners
            w1cen = w1.getCentre();
            w1cor1 = w1.getCorner();
            w1cor2 = -w1cor1;

            w2cen = w2.getCentre();
            w2cor1 = w2.getCorner();
            w2cor2 = -w2cor1;

            % Test 8 intersection points against two wedge centres, and 
            % test if the 4 wedge corners are lying in the other wedge
            vertices = [intn; w1cor1; w1cor2; w2cor1; w2cor2];
            in_wedges = [intn * w1cen' >= 0 & intn * w2cen' >= 0;
                         w1cor1 * w2.n1' >= 0 && w1cor1 * w2.n2' >= 0;
                         w1cor2 * w2.n1' >= 0 && w1cor2 * w2.n2' >= 0;
                         w2cor1 * w1.n1' >= 0 && w2cor1 * w1.n2' >= 0;
                         w2cor2 * w1.n1' >= 0 && w2cor2 * w1.n2' >= 0];

            int_points = vertices(in_wedges, :);
        end
        
        function t_list = wedges2intTPatches(w1, w2, min_size)
            int_points = wedge.wedges2int(w1, w2);
            if size(int_points,1) > 2
                % Set w1 to thinner wedge
                if w1.getWidth() > w2.getWidth() 
                    w = w1; 
                    w1 = w2;
                    w2 = w;
                end
                t_list = w1.wedge2tPatches(min_size, false);

                % Discard patches not lying in other wedge
                satisfy_normal = StereoInterface.angles(t_list.centre_xyz, [w2.n1; w2.n2]) < pi/2 + t_list.thres;
                t_list = t_list.removePatches( ~satisfy_normal(:,1) | ~satisfy_normal(:,2) );
                
            else % No intersection, no patches to search
                t_list = [];
            end
            
        end
        
    end
end

