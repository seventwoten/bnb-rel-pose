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

    [centres, sigma] = subdivide([0,0,0], pi);
    e_max = 0;
    
    queue = cell(1,max_matches); % Set up one block queue per no. of matches 
    for c = 1: size(centres, 1)
        queue{max_matches} = [queue{max_matches} block(centres(c, :), sigma, max_matches)]; 
    end

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
        cen = queue{i}(1).centre; 
        sigma = queue{i}(1).sigma;
        queue{i}(1) = [];
        
        % Given rotation represented by block centre, find angles(Rp, q) 
        R  = aa2mat(cen);
        Rp = (R * p')';

        a = angles(Rp, q);
        
        fprintf("\nIteration %d   sigma: %d pi e_max = %d\n", iter, sigma/pi, e_max);
        
        % Compute block lower bound with stopping threshold
        edges_tight = a < thres_stop;
        C(2:Np+1, 1+Np+1 : 1+Np+Nq) = edges_tight;
        G = digraph(C);
        m_lower = maxflow(G, 1, Np+Nq+2);
        fprintf("Lower bound: %d\n", m_lower);
        
        % Compute block upper bound at sqrt(3)-sigma threshold
        thres = sqrt(3) * sigma;
        if thres > thres_stop
            edges = a < thres;
            C(2:Np+1, 1+Np+1 : 1+Np+Nq) = edges;
            G = digraph(C);
            m_upper = maxflow(G, 1, Np+Nq+2);
            fprintf("Upper bound: %d \t", m_upper);
        else
            m_upper = m_lower;
        end
        
        % Update e_max
        if m_lower > e_max
            fprintf("Replaced e_max %d with %d.\n", e_max, m_lower);
            e_max = m_lower;
            
            % After updating e_max, check queue for discardable blocks
            % (parentUB < e_max)
            fprintf("Discarded blocks below %d!\n", e_max);
            queue([1:e_max-1]) = {[]};
            
        end
        
        % Discard this block, or subdivide/terminate
        if m_upper > 0 && m_upper >= e_max
            if thres > thres_stop
                fprintf("Continue!\n");
                [centres, sigma] = subdivide(cen, sigma);
                for c = 1: size(centres, 1)
                    queue{m_upper} = [queue{m_upper} block(centres(c, :), sigma, m_upper)];
                end
            else
                fprintf("Solution at stopping resolution: [%d %d %d], score: %d-%d, sigma: %d\n", cen, m_lower, m_upper, sigma);
                solutions = [solutions solution(cen, sigma, m_lower, edges_tight)];
            end
        else
            fprintf("Discard!\n");
        end
        
        i = find(~cellfun(@isempty,queue), 1, 'last');
        fprintf("next i: %d\n", i);
        iter = iter + 1;
        
    end
    
    % Filter out solutions below most recent e_max
    for s = size(solutions, 2):-1:1
        if solutions(s).score < e_max
            solutions(s) = [];
        end
    end
    
end

