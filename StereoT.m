classdef StereoT < StereoInterface
    %STEREOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n1
        n2
        t_list          % List of tTriangles
        t_half_len_stop % Stop at this patch size (affects solution accuracy)
        thres_stop_t    % Angle of lower bound cone, from stopping patch size
        
        check_epipole  % option to turn on near-epipole check
        max_edges_p    % reject p points near epipole with edges > max_edges_p
        max_edges_q    % reject q points near epipole with edges > max_edges_q
        
        angleMat1
        angleMat2
        possibleMatches % Np x Nq mask for including only possible matches
        R_e_max         % Current e_max of R outer loop
    end
    
    methods
        function obj = StereoT(p, q, n1, n2, t_list, t_half_len_stop, epipole_threshold, possible_matches, R_e_max)
            %STEREOT Construct an instance of this class
            %   (pass in [] to use default values)
            obj = obj@StereoInterface(p, q);
            obj.n1 = n1;
            obj.n2 = n2;
            obj.t_half_len_stop = t_half_len_stop;
            obj.thres_stop_t = acos(cos(t_half_len_stop)^2);
            
            if isempty(t_list)
                % default t range to search
                obj.t_list = [tPatch([0,pi/2], pi/2), tPatch([pi,pi/2], pi/2)];
            else
                obj.t_list = t_list;
            end
            
            if isempty(epipole_threshold)
                % default: reject those matching >80% of points
                epipole_threshold = 0.8;
            end
            
            if epipole_threshold == -1
                obj.check_epipole = false;
            else
                obj.check_epipole = true;
                obj.max_edges_p = epipole_threshold * obj.Nq;
                obj.max_edges_q = epipole_threshold * obj.Np;
            end
            
            if exist('possible_matches', 'var') && ~isempty(possible_matches)
                assert(size(possible_matches, 1) == obj.Np && size(possible_matches, 2) == obj.Nq, 'PossibleMatches is the wrong size');
                obj.possibleMatches = possible_matches;
            else
                obj.possibleMatches = true(obj.Np, obj.Nq);
            end
            
            if exist('R_e_max', 'var') && ~isempty(R_e_max)
                obj.R_e_max = R_e_max;
            end
        end
        
        function [obj] = setContext(obj, block)
            obj.angleMat1 = obj.findAnglesToNormals(block.centre, obj.n1);
            obj.angleMat2 = obj.findAnglesToNormals(block.centre, obj.n2);
        end
        
        function [block] = updateLowerBound(obj, block) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            positiveRange = pi/2 + obj.thres_stop_t;
            block.edges_LB = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
            
            % Near-epipole check: Remove rows/columns that match too many points
            if obj.check_epipole
                rows = sum(block.edges_LB,2) >= obj.max_edges_p;
                cols = sum(block.edges_LB,1) >= obj.max_edges_q;
                for i = 1:size(block.edges_LB, 3)
                    block.edges_LB(rows(:,:,i), :, i)= 0;
                    block.edges_LB(:, cols(:,:,i), i)= 0;
                end
            end
            block.LB = obj.getMaxBipartiteMatching(block.edges_LB);
        end
        
        function [block] = updateUpperBound(obj, block)
            %UPDATEUPPERBOUND Update block upper bound at threshold(s) in block.thres
            
            % Check if triangle centre is in each wedge
            block.edges_UB = block.edges_LB;
            mask = ~block.edges_LB;
            
            % For wedges that don't contain the centre, check for intersections
            tri_n = [block.n12; block.n23; block.n31];
            intn_n1 = zeros(3, obj.Np * obj.Nq, 3);
            intn_n2 = zeros(3, obj.Np * obj.Nq, 3);
            
            intn_n1(:, mask(:), :) = StereoRT.pairwiseCross(tri_n, reshape(obj.n1(mask(:,:,[1,1,1])), [], 3));
            intn_n1 = reshape(intn_n1, 3, obj.Np, obj.Nq, 3);
            
            intn_n2(:, mask(:), :) = StereoRT.pairwiseCross(tri_n, reshape(obj.n2(mask(:,:,[1,1,1])), [], 3));
            intn_n2 = reshape(intn_n2, 3, obj.Np, obj.Nq, 3);
            
            intn_n1 = permute(intn_n1, [1,4,2,3]);
            intn_n2 = permute(intn_n2, [1,4,2,3]);
            
            intns = [intn_n1; -intn_n1; intn_n2; -intn_n2];
            % From 12 intersections per wedge, find if any are inside both
            % triangle and wedge - check against wedge centre and 3 wedge
            % centres formed by pairs of triangle normals
            
            tri_wedge_cen = [tri_n(1,:) + tri_n(2,:);
                             tri_n(2,:) + tri_n(3,:);
                             tri_n(3,:) + tri_n(1,:)];
            tri_wedge_cen = tri_wedge_cen ./ sqrt(sum(tri_wedge_cen.^2, 2));
            [row, col] = find(mask);
            for i = 1:size(row, 1)
                wedge_cen = squeeze(obj.n1(row(i),col(i), :))' ...
                          + squeeze(obj.n2(row(i),col(i), :))';
                wedge_cen = wedge_cen ./ sqrt(sum(wedge_cen.^2));
                a = obj.angles(intns(:,:,row(i),col(i)), [wedge_cen; tri_wedge_cen]);
                intersect = all(a <= pi/2, 2);
                block.edges_UB(row(i), col(i)) = sum(intersect) >= 2;
            end
            
            block.UB = obj.getMaxBipartiteMatching(block.edges_UB);
        end
        
        function [obj, solutions] = findSolutions(obj, early_stop, parallel_mode)
            % Non-parallel mode
            init_list = obj.t_list;
            if numel(init_list) == 1
                init_list = init_list.subdivide();
            end
            obj = obj.bnb(init_list, obj.thres_stop_t, early_stop); 
            
            solutions = obj.solutions;
        end
        
        function [angleMat] = findAnglesToNormals(obj, t, n)
            %FINDANGLES Find angles between t-vector at patch centre, and normals in n
            [s1,s2,s3] = size(n);
            Nt = size(t, 1);
            n = reshape(n, [], 3);
            
            % Return one angle matrix for each t vector (stacked in 3rd dimension)
            angleMat = reshape(obj.angles(t, n)', s1, s2, Nt);
            angleMat(isnan(angleMat)) = 0;
        end
    end
end

