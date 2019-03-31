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
    
    queue = [];
    for c = 1: size(centres, 1)
        queue = [queue block(centres(c, :), sigma, max_matches)];
    end

    i = 1;
    solutions = [];

    % Prepare capacity matrix
    n_C = 1 + Np + Nq + 1;
    C = zeros(n_C, n_C);
    C(1, 2:Np+1) = 1; % set edges from s to p nodes
    C(1+Np+1 : 1+Np+Nq, end) = 1; % set edges from q nodes to t
    %imagesc(C)
    
    
    while i <= size(queue, 2)
        cen = queue(1,i).centre;
        sigma = queue(1,i).sigma;

        R  = aa2mat(cen);
        Rp = (R * p')';

        a = angles(Rp, q);
        
        fprintf("\nIteration %d   sigma: %d pi e_max = %d\n", i, sigma/pi, e_max);
        
        % Compute block upper bound at sqrt(3)-sigma threshold
        thres = sqrt(3) * sigma;
        edges = a < thres;
        C(2:Np+1, 1+Np+1 : 1+Np+Nq) = edges;
        G = digraph(C);
        m_upper = maxflow(G, 1, Np+Nq+2);
        fprintf("Upper bound: %d \t", m_upper);
        
        % Compute block lower bound with stopping threshold
        edges_tight = a < thres_stop;
        C(2:Np+1, 1+Np+1 : 1+Np+Nq) = edges_tight;
        G = digraph(C);
        m_lower = maxflow(G, 1, Np+Nq+2);
        fprintf("Lower bound: %d\n", m_lower);
        
        % Update e_max
        if m_lower > e_max
            fprintf("Replaced e_max %d with %d.\n", e_max, m_lower);
            e_max = m_lower;
            
            % After updating e_max, check queue for discardable blocks
            % (parentUB < e_max)
            for j = size(queue,2) : -1 : i+1
                if queue(j).parentUB < e_max
                    queue(j) = [];
                    fprintf("Discarded %d!\n", j);
                end
            end
        end
        
        % Discard this block, or subdivide/terminate
        if m_upper > 0 && m_upper >= e_max
            if thres/2 > thres_stop
                fprintf("Continue!\n");
                [centres, sigma] = subdivide(cen, sigma);
                for c = 1: size(centres, 1)
                    queue = [queue block(centres(c, :), sigma, m_upper)];
                end
            else
                fprintf("Solution at stopping resolution: [%d %d %d]\n", cen);
                solutions = [solutions solution(cen, sigma, m_lower, edges)];
            end
        else
            fprintf("Discard!\n");
        end
        
        i = i + 1;
        
    end
    
    % Plot solution matrices
    for s = 1:size(solutions, 2)
        figure, imagesc(solutions(s).edges);
    end
end

