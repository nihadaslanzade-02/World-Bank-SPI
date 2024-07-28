function analyze_spi_data(X)
    K_values = 4:7; % Different values of K
    best_explained_variance = -Inf;
    best_results = struct();
    
    for K = K_values
        % Control dimension
        [n, m] = size(X);
        if K > m
            K = m;
        end
        
        % REDKM
        [Urkm, Arkm, Yrkm, frkm, ~] = REDKM(X, K, 2, 30);
        
        % FKM
        [Ufkm, ~, Yfkm, ffkm, ~] = FKM(X, K, 2, 30);

        %CDPCA
        [~,~,~, Ycdpca,~,~]=CDPCA(X, K, 2, 30);
        
        % Other calculations
        confusion_matrix = Urkm' * Ufkm;
        var_yrkm = sum(var(Yrkm, 1)) / size(Yrkm, 2); % 2 yerine Yrkm'nin kolon sayısı
        var_yfkm = sum(var(Yfkm, 1)) / size(Yfkm, 2); % 2 yerine Yfkm'nin kolon sayısı
        var_ycdpca = sum(var(Ycdpca, 1)) / size(Ycdpca, 2); % 2 yerine Ycdpca'nin kolon sayısı (eğer varsa)
        
        % Rotating factors
        Arkmr = rotatefactors(Arkm);
        
        % K value with the best performance
        explained_variance = max(frkm, ffkm); % Comparing by using explained variance
        if explained_variance > best_explained_variance
            best_explained_variance = explained_variance;
            best_results.K = K;
            best_results.Urkm = Urkm;
            best_results.Ufkm = Ufkm;
            best_results.confusion_matrix = confusion_matrix;
            best_results.var_yrkm = var_yrkm;
            best_results.var_yfkm = var_yfkm;
            best_results.var_ycdpca = var_ycdpca;
            best_results.Arkmr = Arkmr;
        end
    end
    
    % Showing results
    disp('En iyi performans gösteren K değeri:');
    disp(best_results.K);
    disp('Urkm:');
    disp(best_results.Urkm);
    disp('Ufkm:');
    disp(best_results.Ufkm);
    disp('Confusion Matrix:');
    disp(best_results.confusion_matrix);
    disp('var(Yrkm):');
    disp(best_results.var_yrkm);
    disp('var(Yfkm):');
    disp(best_results.var_yfkm);
    disp('var(Ycdpca):');
    disp(best_results.var_ycdpca);
    disp('Arkmr:');
    disp(best_results.Arkmr);
end
