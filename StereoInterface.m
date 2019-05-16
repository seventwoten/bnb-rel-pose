classdef StereoInterface
    %STEREOINTERFACE Abstract class for describing Stereo problems
    
    properties
        p
        q
        Np
        Nq
        
        capacityMat % capacity matrix for use in maxflow
        solutions   % list of solutions
        e_max       % maximum solution score (running max of lower bound)
    end
    
    methods (Abstract)
      obj = setContext(obj, block)
      % Given a bnb block, sets context before computing its bounds
      
      block = updateLowerBound(obj, block)
      % Returns block with new lower bound (todo: remove threshold arg from subclass method)
      
      block = updateUpperBound(obj,block)
      % Returns block with new upper bound (todo: remove threshold arg from subclass method)
      
      obj = findSolutions(obj)
      
    end
    
    methods
        
        function obj = StereoInterface(p, q)
            %BLOCK Construct an instance of this class
            if nargin ~= 0
                obj.p = p;
                obj.q = q;
                obj.Np = size(p, 1);
                obj.Nq = size(q, 1);
            end
        end
        
        function [obj] = bnb(obj, initList, thres_stop, early_stop)
            %BNB Do branch and bound
            
            % Set up variables
            max_matches = min(obj.Np, obj.Nq);

            queue = cell(1,max_matches); % Set up one block queue per no. of matches 
            queue{max_matches} = initList;

            obj.e_max = 0;
            i = max_matches;
            obj.solutions = [];

            iter = 1;
            
            while ~isempty(i)
                add_solution = false;
                
                % dequeue block from front of queue with largest max_matches 
                blk = queue{i}(1);
                queue{i}(1) = [];

                % Do any necessary computation before checking bounds
                obj = obj.setContext(blk);

                fprintf("Iteration %d   sigma: %d pi e_max = %d\n", iter, blk.sigma/pi, obj.e_max);

                % Compute block lower bound with stopping threshold
                blk = obj.updateLowerBound(blk, thres_stop);
                fprintf("Lower bound: %d\n", blk.LB);

                % Compute block upper bound at sqrt(3)-sigma threshold
                if blk.LB == max_matches
                    blk.UB = blk.LB;
                else
                    blk = obj.updateUpperBound(blk, thres_stop);
                end
                fprintf("Upper bound: %d \t", blk.UB);

                % Update e_max
                if blk.LB > obj.e_max
                    fprintf("Replaced e_max %d with %d.\n", obj.e_max, blk.LB);
                    obj.e_max = blk.LB;
                    
                    % Add block to solutions (quick fix for runs that
                    % terminate with no solutions)
                    add_solution = true;
                    
                    % After updating e_max, check queue for discardable blocks
                    % (parentUB < e_max)
                    fprintf("Discarded blocks below %d!\n", obj.e_max);
                    queue([1:obj.e_max-1]) = {[]};

                end

                % Discard this block, or subdivide/terminate
                if blk.UB > 0 && blk.UB >= obj.e_max
                    if blk.thres > thres_stop
                        fprintf("Continue!\n");
                        queue{blk.UB} = [queue{blk.UB} blk.subdivide()];
                    else
                        add_solution = true;
                    end
                else
                    fprintf("Discard!\n");
                end
                
                % Add solution if flag is true
                if add_solution
                    fprintf(['Solution at stopping resolution: [', repmat('%d ',[1,size(blk.centre,2)]), '], score: %d-%d, sigma: %d\n'], blk.centre, blk.LB, blk.UB, blk.sigma);
                    obj.solutions = [obj.solutions blk];
                end
                
                % early_stop option: Terminate early with this solution if e_max reaches max_matches
                if early_stop && obj.e_max == max_matches
                    break;
                end
                
                i = find(~cellfun(@isempty,queue), 1, 'last');
                fprintf("next i: %d\n\n", i);
                iter = iter + 1;

            end

            % Filter out solutions below most recent e_max
            for s = size(obj.solutions, 2):-1:1
                if obj.solutions(s).LB < obj.e_max
                    obj.solutions(s) = [];
                end
            end
            
        end
        
        function [out] = getMaxBipartiteMatching(obj, edges)
            %GETMAXBIPARTITEMATCHING Return size of max bipartite matching
            obj = obj.setCapacityMatrix(obj.Np, obj.Nq, edges);
            
            % Call maxflow to get max bipartite matching
            G = digraph(obj.capacityMat);
            out = maxflow(G, 1, obj.Np + obj.Nq + 2);
        end
        
        function [obj] = setCapacityMatrix(obj, Np, Nq, edges)
            %SETCAPACITYMATRIX Set capacityMat attribute to template matrix
            C = obj.getTemplateCapacityMatrix(Np, Nq);
            C(2:Np+1, 1+Np+1 : 1+Np+Nq) = edges;
            obj.capacityMat = C;
        end
    end
    
    methods (Static)
        function [out] = getTemplateCapacityMatrix(Np, Nq)
            % Prepare shared template capacity matrix
            persistent C;
            if isempty(C)
                n_C = 1 + Np + Nq + 1;
                C = zeros(n_C, n_C);
                C(1, 2:Np+1) = 1; % set edges from s to p nodes
                C(1+Np+1 : 1+Np+Nq, end) = 1; % set edges from q nodes to t
            end
            out = C;
        end
        
        function [result] = angles(v1, v2)
            %ANGLES Finds pairwise angles between two arrays of 3D vectors
            %   Detailed explanation goes here

            v1_mag = sqrt(sum(v1.^2,2)); % N x 1
            v2_mag = sqrt(sum(v2.^2,2)); % N x 1
            divisor = v1_mag * v2_mag';   % N x N

            cos_angle = (v1 * v2') ./ divisor;
            result = abs(acos(cos_angle));
        end

    end
end

