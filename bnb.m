function [solutions] = bnb(p, q, thres_stop)
%BNB Finds camera rotation given two sets of image points 
%   p and q are N x 3 arrays (image points as 3D vectors).
%   thres_stop is the stopping threshold in radians.

%   (Rp-q angles below threshold are considered potential matches. 
%   The threshold shrinks until thres_stop as Bnb progresses.)

    % Set up variables
    Np = size(p, 1);
    Nq = size(q, 1);
    max_matches = min(Np, Nq);

    queue = cell(1,max_matches); % Set up one block queue per no. of matches 
    init_blk = block([0,0,0], pi, max_matches);
    queue{max_matches} = init_blk.subdivide(); 

    e_max = 0;
    i = max_matches;
    solutions = [];

    % Prepare capacity matrix
    n_C = 1 + Np + Nq + 1;
    C = zeros(n_C, n_C);
    C(1, 2:Np+1) = 1; % set edges from s to p nodes
    C(1+Np+1 : 1+Np+Nq, end) = 1; % set edges from q nodes to t
    %imagesc(C)
    
    iter = 1;
    while ~isempty(i)
        % dequeue block from front of queue with largest max_matches 
        blk = queue{i}(1);
        queue{i}(1) = [];
        
        % Given rotation represented by block centre, find angles(Rp, q) 
        R  = aa2mat(blk.centre);
        Rp = (R * p')';

        a = angles(Rp, q);
        
        fprintf("Iteration %d   sigma: %d pi e_max = %d\n", iter, blk.sigma/pi, e_max);
        
        % Compute block lower bound with stopping threshold
        edges_tight = a < thres_stop;
        C(2:Np+1, 1+Np+1 : 1+Np+Nq) = edges_tight;
        G = digraph(C);
        blk.LB = maxflow(G, 1, Np+Nq+2);
        fprintf("Lower bound: %d\n", blk.LB);
        
        % Compute block upper bound at sqrt(3)-sigma threshold
        thres = sqrt(3) * blk.sigma;
        if thres > thres_stop
            edges = a < thres;
            C(2:Np+1, 1+Np+1 : 1+Np+Nq) = edges;
            G = digraph(C);
            blk.UB = maxflow(G, 1, Np+Nq+2);
            fprintf("Upper bound: %d \t", blk.UB);
        else
            blk.UB = blk.LB;
        end
        
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
            if thres > thres_stop
                fprintf("Continue!\n");
                queue{blk.UB} = [queue{blk.UB} blk.subdivide()];
            else
                fprintf("Solution at stopping resolution: [%d %d %d], score: %d-%d, sigma: %d\n", blk.centre, blk.LB, blk.UB, blk.sigma);
                solutions = [solutions solution(blk.centre, blk.sigma, blk.LB, edges_tight)];
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

