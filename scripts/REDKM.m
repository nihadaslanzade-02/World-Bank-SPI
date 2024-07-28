function [Urkm, Arkm, Yrkm, frkm, inrkm] = REDKM(X, K, Q, Rndstart)

    % initialization
    maxiter = 100;
    eps = 1e-10;
    opts.disp = 0;
    VC = eye(Q);

    [n, J] = size(X);

    % centring matrix
    Jm = speye(n) - (1/n) * (sparse(1:n, 1:n, 1, n, n));

    % compute var-covar matrix 
    S = (1/n) * X' * Jm * X;

    % Standardize data
    Xs = Jm * X * diag(diag(S))^-0.5;

    st = sum(sum(Xs.^2));

    for loop = 1:Rndstart
        U = randPU(n, K);
        su = sum(U);
        Xmean = diag(1./su) * U' * Xs;
        it = 0;

        % update A
        XX = Xs' * U * diag(1./su) * U' * Xs;
        [A, L] = eig(XX);
        [dL, idL] = sort(diag(L), 'descend');
        L = diag(dL);
        A = A(:, idL);

        % Determine the number of columns in A
        numColsA = size(A, 2);

        % Ensure Q does not exceed the number of columns in A
        Q = min(Q, numColsA);

        % Select the first Q columns of A
        A = A(:, 1:Q);

        Ymean = Xmean * A;
        Y = Xs * A;
        f0 = trace((1./n) * Ymean' * U' * U * Ymean);

        % iteration phase
        fdif = 2 * eps;
        while fdif > eps && it < maxiter
            it = it + 1;

            % given Ymean update U
            U = zeros(n, K);
            for i = 1:n
                mindif = sum((Xs(i,:) - Ymean(1,:) * A').^2);
                posmin = 1;
                for j = 2:K
                    dif = sum((Xs(i,:) - Ymean(j,:) * A').^2);
                    if dif < mindif
                        mindif = dif;
                        posmin = j;
                    end 
                end
                U(i, posmin) = 1;
            end

            su = sum(U);
            while sum(su == 0) > 0
                [m, p1] = min(su);
                [m, p2] = max(su);
                ind = find(U(:, p2));
                ind = ind(1:floor(su(p2) / 2));
                U(ind, p1) = 1;
                U(ind, p2) = 0;
                su = sum(U);
            end 

            % given U compute Xmean
            Xmean = diag(1./su) * U' * Xs;

            % given U and Xmean update A
            XX = Xs' * U * diag(1./su) * U' * Xs;
            [A, L] = eig(XX);
            [dL, idL] = sort(diag(L), 'descend');
            L = diag(dL);
            A = A(:, idL);

            % Determine the number of columns in A
            numColsA = size(A, 2);

            % Ensure Q does not exceed the number of columns in A
            Q = min(Q, numColsA);

            % Select the first Q columns of A
            A = A(:, 1:Q);

            Ymean = Xmean * A;
            Y = Xs * A;
            f = trace((1./n) * Ymean' * U' * U * Ymean);
            fdif = f - f0;

            if fdif > eps 
                f0 = f; A0 = A; 
            else
                break
            end
        end

        fprintf('REDKM: Loop=%g, Explained variance=%g, iter=%g, fdif=%g\n', loop, f, it, fdif);

        if loop == 1
            Urkm = U;
            Arkm = A;
            Yrkm = Xs * Arkm;
            frkm = f;
            looprkm = 1;
            inrkm = it;
            fdifo = fdif;
        end
        if f > frkm
            Urkm = U;
            frkm = f;
            Arkm = A;
            Yrkm = Xs * Arkm;
            looprkm = loop;
            inrkm = it;
            fdifo = fdif;
        end
    end

    % sort components in descend order of variance and rotate factors
    varYrkm = var(Yrkm, 1);
    [c, ic] = sort(varYrkm, 'descend');
    Arkm = Arkm(:, ic);
    Yrkm = Yrkm(:, ic);

    [c, ic] = sort(diag(Urkm' * Urkm), 'descend');
    Urkm = Urkm(:, ic);
    fprintf('REDKM (Final): Percentage Explained variance=%g, looprkm=%g, iter=%g, fdif=%g\n', frkm / st * 100, looprkm, inrkm, fdifo);
end
