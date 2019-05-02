classdef cube < block
    %CUBE 3D block
    properties
        angleMat
        patches
    end
    
    methods
        function obj = cube(c, s, p)
            %CUBE Construct an instance of this class
            if nargin == 0
                c = [];
                s = [];
                p = [];
            end
            obj = obj@block(c, s, p);
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
                subblocks(c) = cube(centre_new(c,:), sigma_new, obj.UB);
            end
            
        end
        
        function [R] = aa2mat(obj)
            %AA2MAT Returns an R matrix, converted from cube centre as a 3d 
            % point in an R-space ball. (centre is 1 x 3)
            
            alpha = sqrt(sum(obj.centre.^2,2));
            u = obj.centre / alpha;

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

