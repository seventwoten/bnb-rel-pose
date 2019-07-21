classdef StereoT < StereoInterface
    %STEREOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n1
        n2
        thres_stop_t
        t_long_lat
        t_half_len
        
        check_epipole  % option to turn on near-epipole check
        max_edges_p    % reject p points near epipole with edges > max_edges_p
        max_edges_q    % reject q points near epipole with edges > max_edges_q
        
        angleMat1
        angleMat2
    end
    
    methods
        function obj = StereoT(p, q, n1, n2, t_long_lat, t_half_len, t_half_len_stop, epipole_threshold)
            %STEREOT Construct an instance of this class
            %   t_long_lat and t_half_len are optional 
            %   (pass in [] to use default values)
            obj = obj@StereoInterface(p, q);
            obj.n1 = n1;
            obj.n2 = n2;
            obj.thres_stop_t = acos(cos(t_half_len_stop)^2);
            
            if isempty(t_long_lat) || isempty(t_half_len)
                % default t range to search
                obj.t_long_lat = [0, pi];
                obj.t_half_len = pi/2;
            else
                obj.t_long_lat = t_long_lat;
                obj.t_half_len = t_half_len;
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
            
            
        end
        
        function [obj] = setContext(obj, block)
            obj.angleMat1 = obj.findAnglesToNormals(block.centre_xyz, obj.n1);
            obj.angleMat2 = obj.findAnglesToNormals(block.centre_xyz, obj.n2);
        end
        
        function [block] = updateLowerBound(obj, block, thres_stop) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            positiveRange = pi/2 + thres_stop;
            block.edges_stop = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
            
            % Near-epipole check: Remove rows/columns that match too many points
            if obj.check_epipole
                rows = sum(block.edges_stop,2) >= obj.max_edges_p;
                cols = sum(block.edges_stop,1) >= obj.max_edges_q;
                for i = 1:size(block.edges_stop, 3)
                    block.edges_stop(rows(:,:,i), :, i)= 0;
                    block.edges_stop(:, cols(:,:,i), i)= 0;
                end
            end
            for i = 1:size(block.edges_stop, 3)
                block.LB(i) = obj.getMaxBipartiteMatching(block.edges_stop(:,:,i));
            end
        end
        
        function [block] = updateUpperBound(obj, block)
            %UPDATEUPPERBOUND Update block upper bound at threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            positiveRange = pi/2 + block.thres;
            edges = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
            
            for i = 1:size(edges, 3)
                block.UB(i) = obj.getMaxBipartiteMatching(edges(:,:,i));
            end
        end
        
        function [obj, solutions] = findSolutions(obj, early_stop, parallel_mode)
            if ~exist('parallel_mode','var')
                % Default mode: process t blocks in parallel
                parallel_mode = true;
            end
            
            if parallel_mode
                init_blk = tPatchList(obj.t_long_lat, obj.t_half_len);
                obj = obj.bnbList(init_blk.subdivide(), obj.thres_stop_t, early_stop); 
            else
                init_blk = tPatch(obj.t_long_lat, obj.t_half_len);
                obj = obj.bnb(init_blk.subdivide(), obj.thres_stop_t, early_stop); 
            end
            
            solutions = obj.solutions;
        end
        
        function [obj] = bnbList(obj, initList, thres_stop, early_stop)
            %BNBLIST Do branch and bound, processing blocks in parallel using tPatchList
            
            % Set up variables
            max_matches = min(obj.Np, obj.Nq);
            blk = initList;

            obj.e_max = 0;
            obj.solutions = [];

            iter = 1;
            
            while blk.Nb ~= 0
                fprintf("Num of blocks: %d\n", blk.Nb);
                add_solution = false(blk.Nb, 1);
                
                % Do any necessary computation before checking bounds
                obj = obj.setContext(blk);

                fprintf("Iteration %d   sigma: %d pi e_max = %d\n", iter, blk.sigma/pi, obj.e_max);

                % Compute all lower bounds with stopping threshold
                blk = obj.updateLowerBound(blk, thres_stop);
                fprintf(['Lower bound: [', repmat('%d ',[1,size(blk.LB,1)]), '] \n'], blk.LB);
                
                % Compute block upper bound at sqrt(3)-sigma threshold
                if blk.thres <= thres_stop % blk.LB == max_matches ||
                    blk.UB = blk.LB;
                else
                    blk = obj.updateUpperBound(blk);
                end
                fprintf(['Upper bound: [', repmat('%d ',[1,size(blk.UB,1)]), '] \t'], blk.UB);

                % Update e_max
                max_LB = max(blk.LB);
                if max_LB > 0 && max_LB >= obj.e_max
                    fprintf("Replaced e_max %d with %d.\n", obj.e_max, max_LB);
                    obj.e_max = max_LB;
                    
                    % Add blocks to solutions (quick fix for runs that terminate with no solutions)
                    add_solution(blk.LB == obj.e_max) = true;
                end
                
                % Discard or subdivide/terminate current blocks
                blk_cache = blk;
                if blk.thres > thres_stop
                    blk = blk.subdivide(blk.UB > 0 & blk.UB >= obj.e_max);
                else
                    fprintf("Stopping threshold reached\n");
                    add_solution(blk.UB > 0 & blk.UB >= obj.e_max) = true;
                    blk.Nb = 0;
                end
                
                % Add solution from blk_cache if flag is true
                indices = find(add_solution);
                sols = tPatch.empty;
                for i = 1:size(indices,1)
                    n = indices(i);
                    sols(i) = tPatch(blk_cache.centre(n,:),blk_cache.sigma);
                    sols(i).LB = blk_cache.LB(n);
                    sols(i).UB = blk_cache.UB(n);
                    sols(i).edges_stop = blk_cache.edges_stop(:,:,n);
                end
                if ~isempty(sols)
                    obj.solutions = [obj.solutions sols];
                end
                
                % early_stop option: Terminate early with this solution if e_max reaches max_matches
                if early_stop && obj.e_max == max_matches
                    break;
                end
                
                iter = iter + 1;

            end

            % Filter out solutions below most recent e_max
            for s = size(obj.solutions, 2):-1:1
                if obj.solutions(s).LB < obj.e_max
                    obj.solutions(s) = [];
                end
            end
            
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

