for case_number = 1:18
    rep = 20;
    if is_demo==true
        case_number=0;
        rep = 4;
    end
    
    load(sprintf('data/Case%d.mat', case_number));
    fprintf('Calculating unit GLS for Case %d...\n', case_number);
    M = size(X_full{1},2);
    p = size(X_full{1}{1},2);
    q = size(Z_full{1}{1},2);
    n = zeros(M,1);
    for j=1:rep
        fprintf('Replicate %d...\n', j);
        X = X_full{j};
        Z = Z_full{j};
        Y = Y_full{j};
        for i=1:M
            n(i) = size(X{i},1);
        end
        W = cell(1,M);
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Way 1
%         big_Z = zeros(sum(n), M*q);
%         long_Z = zeros(sum(n), q);
%         long_X = zeros(sum(n), p);
%         long_Y = zeros(sum(n),1);
%         group_var = zeros(sum(n),1);
%         for i=1:M
%             big_Z(1+sum(n(1:i-1)):sum(n(1:i)), 1+(i-1)*q:i*q) = Z{i};
%             long_Z(1+sum(n(1:i-1)):sum(n(1:i)), :) = Z{i};
%             long_X(1+sum(n(1:i-1)):sum(n(1:i)), :) = X{i};
%             long_Y(1+sum(n(1:i-1)):sum(n(1:i))) = Y{i};
%             group_var(1+sum(n(1:i-1)):sum(n(1:i))) = i;
%         end
%         lme = fitlmematrix([long_X, big_Z], long_Y, long_Z, group_var, 'CovariancePattern', 'Isotropic','FitMethod','REML');
%         [psi, sigma] = covarianceParameters(lme);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Way 2
        psi = {sigma_u2};
        sigma = sigma_e2;

        for i=1:M
            W{i} = (sigma*eye(n(i))+Z{i}*psi{1}*Z{i}')\eye(n(i));
        end
        
        
        theta_U = zeros(M, q);
        beta_U = zeros(M, p);
        Sigma_big = cell(1,M);
        Var_big = cell(1,M);
        for i=1:M
            T = [X{i},Z{i}];
            Var_big{i} = T'*W{i}*T;
            Sigma_big{i} = Var_big{i} \ eye(p+q);
            unit_GLS = Sigma_big{i} * T'*W{i}*Y{i};
            beta_U(i,:) = unit_GLS(1:p);
            theta_U(i,:) = unit_GLS(p+1:end);
        end
        save(sprintf('data/Case%d_Rep%d_unit_GLS.mat', case_number, j), 'theta_U', 'beta_U', 'Sigma_big', 'W', '-v7.3');
    end
    
    if is_demo==true
        return;
    end
end