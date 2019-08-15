classdef StereoRT < StereoInterface
    %STEREORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_known = 0        % Number of known correspondences
        p_known
        q_known
        possibleMatches    % Np x Nq mask for including only possible matches
        
        R_list             % List of RCubes
        thres_stop_R
        
        t_list             % Default list of tPatches
        t_half_len_stop
        
        delta
        epipole_threshold
        
        n1_LB
        n2_LB
        n1_UB
        n2_UB
        angleMat
    end
    
    methods
        function obj = StereoRT(p, q, R_list, R_sigma_stop, t_list, t_half_len_stop, delta, epipole_threshold, corr_indices, possible_matches)
            %STEREORT Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@StereoInterface(p, q);
            obj.thres_stop_R = sqrt(3) * R_sigma_stop;
            obj.t_half_len_stop = t_half_len_stop;
            obj.delta = delta;
            obj.epipole_threshold = epipole_threshold;
            
            if exist('possible_matches', 'var') && ~isempty(possible_matches)
                assert(size(possible_matches, 1) == obj.Np && size(possible_matches, 2) == obj.Nq, 'PossibleMatches is the wrong size');
                obj.possibleMatches = possible_matches;
            end
            
            % Separate out known correspondence
            if ~isempty(corr_indices)
                obj.n_known = size(corr_indices, 1);
                obj.p_known = p(corr_indices(:,1), :);
                obj.q_known = q(corr_indices(:,2), :);
                obj.p(corr_indices(:,1), :) = [];
                obj.q(corr_indices(:,2), :) = [];
                obj.Np = size(obj.p, 1);
                obj.Nq = size(obj.q, 1);
                if ~isempty(obj.possibleMatches)
                    obj.possibleMatches(corr_indices(:,1), :) = [];
                    obj.possibleMatches(:, corr_indices(:,2)) = [];
                end
            end
            
            if isempty(R_list)
                % default R range to search
                obj.R_list = RCube([0, 0, 0], pi);
            else
                obj.R_list = R_list;
            end
            
            if isempty(t_list)
                % default t range to search
                obj.t_list = [tPatch([0,pi/2], pi/2), tPatch([pi,pi/2], pi/2)];
            else
                obj.t_list = t_list;
            end
            
        end
        
        function [obj] = setContext(obj, block)
            %SETCONTEXT Compute lower bound (n1, n2) from (R, p, q, angles, thres, delta)
            [obj.n1_LB, obj.n2_LB] = getWedges(obj, block, obj.p, obj.q, obj.thres_stop_R);
        end
        
        function [normals1, normals2] = getWedges(obj, block, p, q, threshold_R)
            R  = block.aa2mat();
            Rp = (R * p')';
            obj.angleMat = obj.angles(Rp, q);
            
            e_p = obj.delta + threshold_R;
            e_q = obj.delta;
            
            sin_e_p = sin(e_p);
            sin_e_q = sin(e_q);

            sq_terms = sin_e_p .* sin_e_p + sin_e_q .* sin_e_q;
            coeff_cos = 2 .* sin_e_p .* sin_e_q;
            sin2beta_num = sq_terms + (coeff_cos .* cos(obj.angleMat));
            sin2beta = sin2beta_num ./ ((sin(obj.angleMat)) .^ 2); 
            sin_beta = sqrt(sin2beta); % Np x Nq

            % Pairwise combinations of Rp and q
            y_num = permute(sin_e_p * reshape(q', 1,3,[]) + sin_e_q * Rp, [1,3,2]);
            y_den = sin_beta .* sin(obj.angleMat);

            y = y_num ./ y_den;
            y = y ./ sqrt(sum(y.^2, 3));

            x = obj.pairwiseCross(Rp, q);
            x = x ./ sqrt(sum(x.^2, 3));
            
            z = cross(x, y);
            z = z ./ sqrt(sum(z.^2, 3));

            cos_beta = sqrt(1-sin2beta);

            normals1 = sin_beta .* z + cos_beta .* x;
            normals2 = sin_beta .* z - cos_beta .* x;
            
            % Zero out normals where Rp and q circular ranges overlap 
            % (wedge constraint does not hold)
            has_overlap = obj.angleMat < e_p + e_q;
            normals1(has_overlap(:,:,[1,1,1])) = 0;
            normals2(has_overlap(:,:,[1,1,1])) = 0;
            
            % Zero out normals that are complex (Rp, q cannot form wedge)
            for i = 1:size(normals1,1)
                for j = 1:size(normals1,2)
                    if ~isreal(normals1(i, j, :))
                        normals1(i, j, :) = 0;
                    end
                    if ~isreal(normals2(i, j, :))
                        normals2(i, j, :) = 0;
                    end
                end
            end
        end
        
        function t_list = updateTList(obj, block, t_list, thres_R)
            % Default case: Return given t_list

            % Single known correspondence: Find known p-q wedge, use as new t_list
            if obj.n_known ~= 0
                [n1_wedge, n2_wedge] = obj.getWedges(block, obj.p_known, obj.q_known, thres_R);

                if obj.n_known == 1
                    % If wedge normals are valid, replace t_list with wedge patches
                    if wedge.isValid(n1_wedge, n2_wedge)
                        w = wedge(n1_wedge,n2_wedge);
                        t_list = w.wedge2tPatches(obj.t_half_len_stop);
                    end
                elseif obj.n_known == 2
                    % Check if wedges are valid
                    w1_valid = wedge.isValid(n1_wedge(1,1,:), n2_wedge(1,1,:));
                    w2_valid = wedge.isValid(n1_wedge(2,2,:), n2_wedge(2,2,:));

                    if w1_valid && w2_valid
                        w1 = wedge(n1_wedge(1,1,:), n2_wedge(1,1,:));
                        w2 = wedge(n1_wedge(2,2,:), n2_wedge(2,2,:));
                        t_list = wedge.wedges2intTPatches(w1, w2, obj.t_half_len_stop);
                        
                    elseif w1_valid && ~w2_valid % Only w1 is valid
                        w1 = wedge(n1_wedge(1,1,:), n2_wedge(1,1,:));
                        t_list = w1.wedge2tPatches(obj.t_half_len_stop);
                        
                    elseif ~w1_valid && w2_valid % Only w2 is valid
                        w2 = wedge(n1_wedge(2,2,:), n2_wedge(2,2,:));
                        t_list = w2.wedge2tPatches(obj.t_half_len_stop);
                        
                    end
                end
            end
        end
        
        function [block] = updateLowerBound(obj, block) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.n1_LB) & ~isempty(obj.n2_LB), 'Context was not set');
            
            % Update T search list if known correspondence exists
            tlist = obj.updateTList(block, obj.t_list, obj.thres_stop_R);
            
            if ~isempty(tlist)
                % Pass near-epipole check option to T search, only at R stopping threshold
                st = StereoT(obj.p, obj.q, obj.n1_LB, obj.n2_LB, tlist, obj.t_half_len_stop, obj.epipole_threshold, obj.possibleMatches);
                fprintf("{\n");
                st = st.findSolutions(true); % early_stop = true
                fprintf("}\n");
                block.LB = st.e_max;
                block.patches = st.solutions;
            else
                fprintf("Intersection of wedges is empty.\nLower bound: 0\n");
                block.LB = 0;
            end
        end
        
        function [block] = updateUpperBound(obj, block)
            %UPDATEUPPERBOUND Update block upper bound at threshold
            [obj.n1_UB, obj.n2_UB] = getWedges(obj, block, obj.p, obj.q, block.thres);
            
            % Update T search list if known correspondence exists
            tlist = obj.updateTList(block, obj.t_list, block.thres);

            if ~isempty(tlist)
                st = StereoT(obj.p, obj.q, obj.n1_UB, obj.n2_UB, tlist, obj.t_half_len_stop, -1, obj.possibleMatches);
                fprintf("{\n");
                st = st.findSolutions(true); % early_stop = true
                fprintf("}\n");
                block.UB = st.e_max;
            else
                fprintf("Intersection of wedges is empty.\nUpper bound: 0\n");
                block.UB = 0;
            end
        end
        
        function [obj, solutions] = findSolutions(obj, early_stop)
            % Default: Return all surviving blocks as solutions
            if ~exist('early_stop', 'var') || isempty(early_stop)
                early_stop = false;
            end
            
            init_list = obj.R_list;
            if numel(init_list) == 1
                init_list = init_list.subdivide();
            end
            obj = obj.bnb(init_list, obj.thres_stop_R, early_stop);
            solutions = obj.solutions;
        end
        
    end
    
    methods (Static)
        function [out] = pairwiseCross(v1, v2)
        %PAIRWISECROSS Returns pairwise cross products of two arrays of
        %vectors. Arrays are (num_vectors x dimensions)
            n1 = size(v1, 1);
            n2 = size(v2, 1);
            dim = size(v1,2);
            out = zeros(n1, n2, dim);

            for row = 1 : n1
                out(row,:,:) = cross(repmat(v1(row, :), n2, 1), v2);
            end
        end
        
    end
end

