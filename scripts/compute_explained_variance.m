function explainedVariances = compute_explained_variance(X, K_values)
    % This function computes the explained variance for given K values
    % Input:
    %   X - the dataset (each row is a data point)
    %   K_values - a vector containing the K values to compute explained variances for
    % Output:
    %   explainedVariances - a vector of explained variances for each K value

    % Compute the overall mean of the dataset
    overallMean = mean(X);

    % Compute the total sum of squares (TSS)
    TSS = sum(sum(bsxfun(@minus, X, overallMean).^2));

    % Initialize the explained variances array
    explainedVariances = zeros(length(K_values), 1);

    % Loop over each K value
    for idx = 1:length(K_values)
        K = K_values(idx);

        % Perform k-means clustering using the optimized kmeansVICHI function
        [~, U, ~, ~] = kmeansN(X, K, 10);

        % Compute the within-cluster sum of squares (WCSS)
        centroids = bsxfun(@rdivide, U' * X, sum(U)');
        WCSS = sum(arrayfun(@(j) sum(sum(bsxfun(@minus, X(U(:, j) == 1, :), centroids(j, :)).^2)), 1:K));

        % Compute the explained variance
        explainedVariances(idx) = 1 - WCSS/TSS;

         % Display the explained variances for each K
        fprintf('K = %d: %f\n', K, explainedVariances(idx));
    end
    
    % Plot the explained variances
    figure;
    plot(K_values, explainedVariances, '-o');
    xlabel('Number of Clusters (K)');
    ylabel('Explained Variance');
    title('Explained Variance vs. Number of Clusters');
    grid on;
end