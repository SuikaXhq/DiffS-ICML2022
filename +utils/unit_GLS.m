function [theta, psi, sigma] = unit_GLS(X, Z, Y)
    % Calculate unit-wise GLS estimators
    
    % Input:
    % X: 1xM Cell with n_i x p Matrix contents
    % Z: 1xM Cell with n_i x q Matrix contents
    % Y: 1xM Cell with n_i-d Vector contents
    
    % Output:
    % theta: 
    % 

    M = size(X,2);
    p = size(X{1},2);
    q = size(Z{1},2);
    n = zeros(M,1);
    for i=1:M
        n(i) = size(X{i},1);
    end
    %fprintf('Calculating W_i..\n');
    W = cell(1,M);
    big_Z = zeros(sum(n), M*q);
    long_Z = zeros(sum(n), q);
    long_X = zeros(sum(n), p);
    long_Y = zeros(sum(n),1);
    group_var = zeros(sum(n),1);
    for i=1:M
        big_Z(1+sum(n(1:i-1)):sum(n(1:i)), 1+(i-1)*q:i*q) = Z{i};
        long_Z(1+sum(n(1:i-1)):sum(n(1:i)), :) = Z{i};
        long_X(1+sum(n(1:i-1)):sum(n(1:i)), :) = X{i};
        long_Y(1+sum(n(1:i-1)):sum(n(1:i))) = Y{i};
        group_var(1+sum(n(1:i-1)):sum(n(1:i))) = i;
    end
    lme = fitlmematrix([long_X, big_Z], long_Y, long_Z, group_var, 'CovariancePattern', 'Isotropic','FitMethod','REML')
    [psi, sigma] = covarianceParameters(lme);
    theta = 0;
%     for i=1:M
%         W{i} = (sigma*eye(n(i))+Z{i}*psi{1}*Z{i}')\eye(n(i));
%     end
%     %fprintf('Initialization done.\n');
% 
%     %% Step 1: Calculate unit GLS estimates
%     %fprintf('Step 1: Calculate unit GLS estimates.\n');
%     theta_U = zeros(M, q);
%     Sigma_big = cell(1,M);
%     Var_big = cell(1,M);
%     tic;
%     for i=1:M
%         T = [X{i},Z{i}];
%         Var_big{i} = T'*W{i}*T;
%         Sigma_big{i} = Var_big{i} \ eye(p+q);
%         theta_U(i,:) = Sigma_big{i}(p+1:end, :) * T'*W{i}*Y{i};
%     end
%     timecost(1) = toc;
%     %fprintf('Step 1 done. Timecost: %.6fs\n',timecost(1));