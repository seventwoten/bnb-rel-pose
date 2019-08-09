classdef StereoT < StereoInterface
    %STEREOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n1
        n2
        t_list          % List of tPatches
        t_half_len_stop % Stop at this patch size (affects solution accuracy)
        thres_stop_t    % Angle of lower bound cone, from stopping patch size
        
        check_epipole  % option to turn on near-epipole check
        max_edges_p    % reject p points near epipole with edges > max_edges_p
        max_edges_q    % reject q points near epipole with edges > max_edges_q
        
        angleMat1
        angleMat2
        possibleMatches % Np x Nq mask for including only possible matches
    end
    
    methods
        function obj = StereoT(p, q, n1, n2, t_list, t_half_len_stop, epipole_threshold, possible_matches)
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
            end
            
        end
        
        function [obj] = setContext(obj, block)
            obj.angleMat1 = obj.findAnglesToNormals(block.centre_xyz, obj.n1);
            obj.angleMat2 = obj.findAnglesToNormals(block.centre_xyz, obj.n2);
        end
        
        function [block] = updateLowerBound(obj, block) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            positiveRange = pi/2 + obj.thres_stop_t;
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
            
            % Filter to include possible matches only
            if ~isempty(obj.possibleMatches)
                block.edges_stop = block.edges_stop & obj.possibleMatches;
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
            
            % Filter to include possible matches only
            if ~isempty(obj.possibleMatches)
                edges = edges & obj.possibleMatches;
            end
            
            for i = 1:size(edges, 3)
                block.UB(i) = obj.getMaxBipartiteMatching(edges(:,:,i));
            end
        end
        
        function [obj, solutions] = findSolutions(obj, early_stop, parallel_mode)
            if ~exist('parallel_mode','var') || isa(obj.t_list,'tPatchList')
                % Default mode: process t blocks in parallel
                parallel_mode = true;
            end
            
            if parallel_mode
                if isa(obj.t_list,'tPatchList')
                    init_list = obj.t_list;
                else % Wrap list of tPatch objects into tPatchList
                    Nt = numel(obj.t_list);
                    if Nt == 1
                        init_list = tPatchList(obj.t_list.centre, obj.t_list.sigma);
                        init_list = init_list.subdivide();
                    else
                        centres = zeros(Nt, 2);
                        sigmas  = zeros(Nt, 1);
                        for i = 1:Nt
                            % Extract centres and sigmas from list of tPatches
                            centres(i, :) = obj.t_list(i).centre;
                            sigmas(i)     = obj.t_list(i).sigma;
                        end

                        % Enforce equal sigmas
                        assert(all(abs(sigmas - sigmas(1)) < 1e-8), 'Size of initial t-patches must be the same');
                        init_list = tPatchList(centres, sigmas(1));
                    end
                end
                
                obj = obj.bnbList(init_list, obj.thres_stop_t, early_stop); 
            else % Non-parallel mode
                init_list = obj.t_list;
                if numel(init_list) == 1
                    init_list = init_list.subdivide();
                end
                obj = obj.bnb(init_list, obj.thres_stop_t, early_stop); 
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
                blk = obj.updateLowerBound(blk);
                fprintf(['Lower bound: [', repmat('%d ',[1,min(numel(blk.LB), 20)]), '] \n'], blk.LB(1: min(numel(blk.LB), 20))); 
                
                % Compute block upper bound at sqrt(3)-sigma threshold
                if blk.sigma <= obj.t_half_len_stop % blk.LB == max_matches ||
                    blk.UB = blk.LB;
                else
                    blk = obj.updateUpperBound(blk);
                end
                fprintf(['Upper bound: [', repmat('%d ',[1,min(numel(blk.UB), 20)]), '] '], blk.UB(1: min(numel(blk.UB), 20)));

                % Update e_max
                max_LB = max(blk.LB);
                if max_LB > 0 && max_LB >= obj.e_max
                    fprintf("\tReplaced e_max %d with %d.", obj.e_max, max_LB);
                    obj.e_max = max_LB;
                    
                    % Add blocks to solutions (quick fix for runs that terminate with no solutions)
                    add_solution(blk.LB == obj.e_max) = true;
                end
                
                % Discard or subdivide/terminate current blocks
                blk_cache = blk;
                surviving_blocks =  blk.UB > 0 & blk.UB >= obj.e_max;
                if blk.sigma > obj.t_half_len_stop
                    if any(surviving_blocks)
                        blk = blk.subdivide(surviving_blocks);
                    else % No more subdivision - terminate
                        blk.Nb = 0;
                    end
                else
                    fprintf("\tStopping threshold reached");
                    add_solution(surviving_blocks) = true;
                    blk.Nb = 0;
                end
                fprintf("\n");
                
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

