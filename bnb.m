function [solutions] = bnb(init_blk, p, q, thres_stop)
%BNB Finds camera rotation given two sets of image points 
%   p and q are N x 3 arrays (image points as 3D vectors).
%   thres_stop is the stopping threshold in radians.

%   (Rp-q angles below threshold are considered potential matches. 
%   The threshold shrinks until thres_stop as Bnb progresses.)

    % Set up variables
    max_matches = min(size(p, 1), size(q, 1));

    queue = cell(1,max_matches); % Set up one block queue per no. of matches 
    queue{max_matches} = init_blk.subdivide(); 

    e_max = 0;
    i = max_matches;
    solutions = [];
    
    iter = 1;
    while ~isempty(i)
        % dequeue block from front of queue with largest max_matches 
        blk = queue{i}(1);
        queue{i}(1) = [];
        
        % Given rotation represented by block centre, find angles(Rp, q) 
        blk = blk.setContext(p, q);
        
        fprintf("Iteration %d   sigma: %d pi e_max = %d\n", iter, blk.sigma/pi, e_max);
        
        % Compute block lower bound with stopping threshold
        blk = blk.setLowerBound(thres_stop);
        fprintf("Lower bound: %d\n", blk.LB);
        
        % Compute block upper bound at sqrt(3)-sigma threshold
        blk = blk.setUpperBound(thres_stop);
        fprintf("Upper bound: %d \t", blk.UB);
        
        % Update e_max
        if blk.LB > e_max
            fprintf("Replaced e_max %d with %d.\n", e_max, blk.LB);
            e_max = blk.LB;
            
            % After updating e_max, check queue for discardable blocks
            % (parentUB < e_max)
            fprintf("Discarded blocks below %d!\n", e_max);
            queue([1:e_max-1]) = {[]};
            
        end
        
        % Discard this block, or subdivide/terminate
        if blk.UB > 0 && blk.UB >= e_max
            if blk.thres > thres_stop
                fprintf("Continue!\n");
                queue{blk.UB} = [queue{blk.UB} blk.subdivide()];
            else
                fprintf("Solution at stopping resolution: [%d %d %d], score: %d-%d, sigma: %d\n", blk.centre, blk.LB, blk.UB, blk.sigma);
                solutions = [solutions solution(blk.centre, blk.sigma, blk.LB, blk.edges_tight)];
            end
        else
            fprintf("Discard!\n");
        end
        
        i = find(~cellfun(@isempty,queue), 1, 'last');
        fprintf("next i: %d\n\n", i);
        iter = iter + 1;
        
    end
    
    % Filter out solutions below most recent e_max
    for s = size(solutions, 2):-1:1
        if solutions(s).score < e_max
            solutions(s) = [];
        end
    end
    
end

