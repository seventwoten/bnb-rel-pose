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
        
        function [obj] = setContext(obj, context)
            %SETCONTEXT Necessary computation before checking bounds
            % Given rotation at block centre, find angles(Rp, q)
            R  = obj.aa2mat();
            Rp = (R * context.p')';
            obj.angleMat = angles(Rp, context.q);
        end
        
        function [obj] = setLowerBound(obj, thres_stop)
            %SETLOWERBOUND Set lower bound at stopping threshold
            assert(~isempty(obj.angleMat), 'Context was not set');
            obj.edges_stop = obj.angleMat < thres_stop;
            obj.LB = obj.getMaxBipartiteMatching(obj.edges_stop);
        end
        
        function [obj] = setUpperBound(obj, thres_stop)
            %SETUPPERBOUND Set upper bound at sqrt(3) * sigma threshold
            assert(~isempty(obj.angleMat), 'Context was not set');
            if obj.thres > thres_stop
                edges = obj.angleMat < obj.thres;
                obj.UB = obj.getMaxBipartiteMatching(edges);
            else
                obj.UB = obj.LB;
            end
        end
        
    end
    
end

