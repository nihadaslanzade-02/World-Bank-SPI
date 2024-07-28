% Defining the range of K values to evaluate
K_values = 2:6;  
pF_values = zeros(size(K_values));  % Preallocating array for pF values

for i = 1:length(K_values)
    K = K_values(i);
    % Running K-means with REDKM function
    [Urkm, Arkm, Yrkm, frkm, inrkm] = REDKM(X, K, 4, 100);
    
    % Calculating pseudo F-statistic using psF function
    [pF, Dw, Db] = psF(X, Urkm);
    
    % Storing the pF value
    pF_values(i) = pF;
end

% Plotting the pF values against K
figure;
plot(K_values, pF_values, 'bo-'); 
xlabel('Number of Clusters (K)');
ylabel('Pseudo F-statistic (pF)');
title('Elbow Method for Determining Optimal K');
grid on;