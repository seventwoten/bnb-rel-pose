classdef RCube < block
    %RCUBE 3D block
    properties
        angleMat
        patches
    end
    
    methods
        function obj = RCube(c, s)
            %RCUBE Construct an instance of this class
            if nargin == 0
                c = [];
                s = [];
            end
            obj = obj@block(c, s);
            obj.thres = sqrt(3) * obj.sigma;
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
                subblocks(c) = RCube(centre_new(c,:), sigma_new);
            end
            
        end
        
        function [R] = aa2mat(obj)
            %AA2MAT Returns an R matrix, converted from cube centre as a 3d 
            % point in an R-space ball. (centre is 1 x 3)
            if sum(abs(obj.centre)) < 1e-06
                R = eye(3);
            else
                alpha = sqrt(sum(obj.centre.^2,2));
                u = obj.centre ./ alpha;

                c = cos(alpha);
                s = sin(alpha);
                t = 1-c;

                U_outer = u' * u;
                U_as = [      0, -u(1,3),  u(1,2) ; 
                         u(1,3),       0, -u(1,1) ; 
                        -u(1,2),  u(1,1),      0 ];

                R = c * eye(3) + t * U_outer + s * U_as;
            end
        end
        
    end
    methods (Static)
        function [R] = rpy2R(roll, pitch, yaw)
            %RPY2R Converts roll, pitch, yaw (right-handed and in radians) to a rotation matrix.

            cos_r = cos(roll);
            cos_p = cos(pitch);
            cos_y = cos(yaw);

            sin_r = sin(roll);
            sin_p = sin(pitch);
            sin_y = sin(yaw);

            R = [cos_r*cos_y,  cos_r*sin_y*sin_p - sin_r*cos_p,  cos_r*sin_y*cos_p + sin_r*sin_p ;
                 sin_r*cos_y,  sin_r*sin_y*sin_p + cos_r*cos_p,  sin_r*sin_y*cos_p - cos_r*sin_p ;
                 -sin_y     ,  cos_y*sin_p                    ,  cos_y*cos_p                    ];
        end
        
        function [rpy] = R2rpy(R)
            %R2RPY Converts rotation matrix to roll, pitch, yaw (right-handed and in radians).

            sy = sqrt(R(1,1) * R(1,1) + R(2,1) * R(2,1));
            singular = sy < 1e-6;

            if ~singular
                r = atan2( R(2,1), R(1,1));
                p = atan2( R(3,2), R(3,3));
                y = atan2(-R(3,1), sy);
            else
                p = atan2(-R(2,3), R(2,2));
                y = atan2(-R(3,1), sy);
                r = 0;
            end

            rpy = [r, p, y];
        end
        
        function [u] = R2aa(R)
            %R2AA Converts an R matrix to a scaled axis-angle
            %representation (vector direction is axis, magnitude is angle) 
            if sum(sum(abs(R - eye(3)))) < 1e-6
                u = [0,0,0];
            else
                u3sin_a = (R(2,1) - R(1,2)) / 2;
                u1sin_a = (R(3,2) - R(2,3)) / 2;
                u2sin_a = (R(1,3) - R(3,1)) / 2;
                
                u = [u1sin_a, u2sin_a, u3sin_a]; 
                abs_sin_a = sqrt(sum(u.^2));
                u = u ./ abs_sin_a;
                
                % u direction is such that alpha is between -pi and pi
                if u(1) > 1e-6 || u(2) > 1e-6
                    cos_a = 1 - (R(1,1) - R(2,2))/(u(1)*u(1) - u(2)*u(2));
                else % u(3) > 1e-6
                    cos_a = 1 - (R(1,1) - R(3,3))/(u(1)*u(1) - u(3)*u(3));
                end
                
                alpha = atan2(abs_sin_a, cos_a);
                if alpha <= -pi || alpha > pi
                    u = -u;
                    alpha = atan2(-abs_sin_a, cos_a);
                end
                
                u = u .* alpha;
            end
        end
        
        function [R_align] = align2R(vec1, vec2)
            % ALIGN2R Find rotation that aligns vec1 to vec2 (both 1x3)
            v   = cross(vec1, vec2);
            v_x = [   0, -v(3),  v(2); 
                   v(3),     0, -v(1); 
                  -v(2),  v(1),     0];

            c = vec1 * vec2';
            if abs(c-1) < 1e-8 || abs(c+1) < 1e-8 % No rotation needed
                R_align = eye(3);
            else
                R_align = eye(3) + v_x + (v_x * v_x ./ (1+c));
            end
        end

    end
end

