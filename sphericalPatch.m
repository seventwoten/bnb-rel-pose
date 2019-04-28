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
            obj.thres = acos(cos(s)^2);
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
        
        function [obj] = setContext(obj, context)
            obj.angleMat1 = obj.angles(context.n1);
            obj.angleMat2 = obj.angles(context.n2);
        end
        
        function [angleMat] = angles(obj, n)
            %ANGLES Find angles between t-vector at patch centre, and normals in n
            [s1,s2,s3] = size(n);
            n = reshape(n, [], 3);
            obj.centre_xyz = obj.spherical2Cartesian(1, obj.centre(1), obj.centre(2));
            angleMat = reshape(angles(obj.centre_xyz, n), s1, s2);
        end
        
        function [obj] = setLowerBound(obj, thres_stop)
            %SETLOWERBOUND Set lower bound at stopping threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            positiveRange = pi/2 + thres_stop;
            obj.edges_stop = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
            obj.LB = obj.getMaxBipartiteMatching(obj.edges_stop);
        end
        
        function [obj] = setUpperBound(obj, thres_stop)
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            if obj.thres > thres_stop
                positiveRange = pi/2 + obj.thres;
                edges = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
                obj.UB = obj.getMaxBipartiteMatching(edges);
            else
                obj.UB = obj.LB;
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


