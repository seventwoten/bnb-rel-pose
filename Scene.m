classdef Scene
    %SCENE Set-up of 3d scene points and two views 
    %   E.g. s = Scene([], 30); s = s.setCam2(); s = s.setView1(); s = s.setView2();
    %   Image vectors are stored in view1 and view2 properties.
    properties
        points3d
        N

        % 3D scene range (length-2 cube, centred at [0,0,2])
        len3d = 2
        cen3d = [0,0,2]
        
        % Camera range (if camera 2 is randomly positioned)
        cam_range_cen   = [0,0,2]
        cam_range_rad   = 2
        cam_range_theta = [0, 2*pi] 
        cam_range_phi   = [0, 2*pi] %[pi/2, 3/4 * pi] % Second camera is 45-90 deg away in phi dimension
        
        % Camera parameters
        f  = 1              % focal length
        bu = 1, bv = 1      % scaling factors
        u0 = 0, v0 = 0      % offsets

        % First camera position (0, 0, 0) and orientation (aligned to world, rows are i, j, k vectors)
        cam1_xyz = [0, 0, 0]
        cam1_R  = eye(3)

        % Second camera position and orientation
        cam2_xyz = [0, 0, 0]    % Cartesian
        cam2_tp  = [0, 0]       % Spherical
        
        cam2_R  = eye(3)
        cam2_rpy = [0, 0, 0]
        cam2_aa  = [0, 0, 0]

        % Image vector noise standard deviation (radians)
        noise_sd = 5.7596e-04
        
        % 2d views of 3d scene
        view1
        view2
    end
    
    methods
        function obj = Scene(points3d, n)
            %SCENE Construct an instance of this class
            %   Either a matrix of 3d points (n x 3) or n must be supplied. 
            if isempty(points3d)
                obj.N = n;
                obj.points3d = rand(n, 3)*obj.len3d - obj.len3d/2 + obj.cen3d;
            else
                obj.points3d = points3d;
                obj.N = size(points3d, 1);
            end
        end
        
        function obj = setCam2(obj, xyz, rpy)
            %SETCAM1 Sets xyz-rpy params of camera 2, relative to camera 1 at (0,0,0,0,0,0)
            if exist('xyz', 'var') && exist('rpy', 'var') && ~isempty(xyz) && ~isempty(rpy) 
                obj.cam2_xyz = obj.cam1_xyz + xyz;
                rtp = tPatch.cartesian2Spherical(obj.cam2_xyz(1),obj.cam2_xyz(2),obj.cam2_xyz(3));
                obj.cam2_tp  = [rtp(2), rtp(3)]; 
                obj.cam2_rpy = rpy;
                obj.cam2_R  = RCube.rpy2R(rpy(1),rpy(2),rpy(3));
                obj.cam2_aa = RCube.R2aa(obj.cam2_R);
            else
                % Position second camera randomly in given ranges, at cam_range_rad from scene centre
                theta_scene = rand() * (obj.cam_range_theta(2) - obj.cam_range_theta(1)) + obj.cam_range_theta(1);
                phi_scene   = rand() * (  obj.cam_range_phi(2) -   obj.cam_range_phi(1)) + obj.cam_range_phi(1);
                obj.cam2_xyz = tPatch.spherical2Cartesian(obj.cam_range_rad, theta_scene, phi_scene) + obj.cam_range_cen;
                rtp          = tPatch.cartesian2Spherical(obj.cam2_xyz(1),obj.cam2_xyz(2),obj.cam2_xyz(3));
                obj.cam2_tp  = [rtp(2), rtp(3)];
                
                % Find rotation to point camera at centre of scene, [0,0,2]
                cam_pos_scene = obj.cam2_xyz - obj.cen3d;
                cam_pos2_unit = cam_pos_scene / sqrt(sum(cam_pos_scene.^2));
                v   = cross([0,0,1], -cam_pos2_unit);
                v_x = [   0, -v(3),  v(2); 
                       v(3),     0, -v(1); 
                      -v(2),  v(1),     0];

                c = [0,0,1] * -cam_pos2_unit';
                if abs(c-1) < 1e-8 || abs(c+1) < 1e-8
                    obj.cam2_R = eye(3);
                else
                    obj.cam2_R = eye(3) + v_x + (v_x * v_x ./ (1+c));
                end
                
                obj.cam2_rpy = RCube.R2rpy(obj.cam2_R);
                obj.cam2_aa  = RCube.R2aa(obj.cam2_R);
            end
        end
        
        function obj = setView1(obj, noise_sd)
            %SETVIEW1 Sets point coords in first camera view (view1),
            % with Gaussian noise in pixel positions of N(0, noise_sd^2). 
            
            [uv] = obj.projectPoints(obj.cam1_xyz, obj.cam1_R);
            obj.view1 = [uv, repmat(obj.f, size(uv, 1), 1)];
            
            % Add Gaussian noise
            if ~exist('noise_sd', 'var') || isempty(noise_sd)
                noise_sd = obj.noise_sd;
            end
            obj.view1 = obj.addNoise(obj.view1, noise_sd);
        end
        
        function obj = setView2(obj, noise_sd)
            %SETVIEW2 Sets point coords in second camera view (view2),
            % with Gaussian noise in pixel positions of N(0, noise_sd^2). 
            
            [uv] = obj.projectPoints(obj.cam2_xyz, obj.cam2_R);
            obj.view2 = [uv, repmat(obj.f, size(uv, 1), 1)];
            
            % Add Gaussian noise
            if ~exist('noise_sd', 'var') || isempty(noise_sd)
                noise_sd = obj.noise_sd;
            end
            obj.view2 = obj.addNoise(obj.view2, noise_sd);
        end
        
        function [uv] = projectPoints(obj, cam_pos, cam_or)
            % PROJECTPOINTS Returns 2D projections of points3d (N x 2), 
            % given camera position and orientation  
            
            % Compute extrinsic matrix to tranform points from world coordinates to camera coordinates
            R = cam_or' ;
            t = -R * cam_pos' ;
            extrinsic = [R t]; 
            world_coords = [obj.points3d ones(obj.N, 1)];
            
            % Scene points in camera coordinates
            s = (extrinsic * world_coords')';
            i_f = [1.0,0.0,0.0]';
            j_f = [0.0,1.0,0.0]';
            k_f = [0.0,0.0,1.0]';

            u = (obj.f * obj.bu * (s * i_f)) ./ (s * k_f) + obj.u0;
            v = (obj.f * obj.bv * (s * j_f)) ./ (s * k_f) + obj.v0;
            
            % Set points behind camera to [inf, inf]
            u(s(:, 3) < 0) = inf;
            v(s(:, 3) < 0) = inf;
            uv = [u v];
        end
    end
    
    methods (Static)
        
        function [img_n] = addNoise(img, noise_sd)
            % Adds angular noise to image points with std. dev. noise_sd
            
            img_n = zeros(size(img));
            for i = 1:size(img,1)
                % For each image vector, find a perpendicular rotation axis
                xy = rand(1,2);
                z = -( xy(1)*img(i,1) + xy(2)*img(i,2) ) ./ img(i,3);
                axis = [xy, z];
                
                % Scale the axis by the magnitude of the noise
                noise = random('Normal', 0, noise_sd);
                axis = axis ./ sqrt(sum(axis.^2)) .* noise;
                
                % Convert to rotation matrix to perturb the vector
                R = RCube(axis, 1).aa2mat();
                img_n(i,:) = R * img(i, :)';
            end
            
        end
        
    end

end

