function [loopOtt,UOtt,fOtt,iterOtt] = kmeansN(X, K, Rndstart)
    % n = number of objects
    % J = number of variables
    % K = number of clusters of the partition
    % maxiter = max number of iterations
    
    maxiter = 100;
    n = size(X,1);
    J = size(X,2);
    epsilon = 0.000001;
    
    for loop = 1:Rndstart
        U0 = randPU(n, K);
        su = sum(U0);
        
        % Given U, compute Xmean (compute centroids)
        Xmean0 = bsxfun(@rdivide, U0' * X, su');
        
        for iter = 1:maxiter
            % Assign each unit to the closest cluster
            U = zeros(n, K);
            for i = 1:n
                distances = sum(bsxfun(@minus, X(i,:), Xmean0).^2, 2);
                [~, posmin] = min(distances);
                U(i, posmin) = 1;
            end
            
            % Given a partition of units, compute Xmean (compute centroids)
            su = sum(U);
            while any(su == 0)
                [~, p1] = min(su);
                [~, p2] = max(su);
                ind = find(U(:,p2));
                ind = ind(1:floor(su(p2)/2));
                U(ind,p1) = 1;
                U(ind,p2) = 0;
                su = sum(U);
            end
            
            % Compute new centroids
            Xmean = bsxfun(@rdivide, U' * X, su');
            
            % Compute objective function
            BB = U * Xmean - X;
            f = sum(sum(BB.^2));
            
            % Stopping rule
            dif = sum(sum((Xmean - Xmean0).^2));
            if dif > epsilon  
                Xmean0 = Xmean;
            else
                break;
            end
        end
        
        if loop == 1
            UOtt = U;
            fOtt = f;
            loopOtt = 1;
            iterOtt = iter;
            XmeanOtt = Xmean;
        end
        fprintf('k-means: loopOtt=%g, fOtt=%g, iter=%g\n', loopOtt, fOtt, iter);
        
        if f < fOtt
            UOtt = U;
            fOtt = f;
            loopOtt = loop;
            iterOtt = iter;
            XmeanOtt = Xmean;
        end
    end
    fprintf('k-means: loopOtt=%g, fOtt=%g, iter=%g\n', loopOtt, fOtt, iter);
end

function [U] = randPU(n, c)
    % Generates a random partition of n objects in c classes
    U = zeros(n, c);
    U(1:c, :) = eye(c);
    U(c+1:n, 1) = 1;
    for i = c+1:n
        U(i, [1:c]) = U(i, randperm(c));
    end
    U(:,:) = U(randperm(n),:);
end
