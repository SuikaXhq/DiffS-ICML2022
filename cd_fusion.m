function [beta, alpha, theta, subgroup, timecost, lambda_full] = cd_fusion(X, Z, Y, train_idx, valid_idx, beta_check, theta_check, W)
%% Initialize
% fprintf('Initializing..\n');
if nargin >= 5
    X_val = X(valid_idx);
    Z_val = Z(valid_idx);
    Y_val = Y(valid_idx);
    X = X(train_idx);
    Z = Z(train_idx);
    Y = Y(train_idx);
    if nargin >= 8
        beta_check_val = beta_check(valid_idx,:);
        theta_check_val = theta_check(valid_idx,:);
        W_val = W(valid_idx);
        beta_check = beta_check(train_idx,:);
        theta_check = theta_check(train_idx,:);
        W = W(train_idx);
    end
end

epsilon = 1e-6;
timecost = zeros(1,6);
M = size(X,2);
p = size(X{1},2);
q = size(Z{1},2);
n = zeros(M,1);
for i=1:M
    n(i) = size(X{i},1);
end
N = sum(n);
% fprintf('M = %d, p = %d, q = %d, N = %d\n', M, p, q, sum(n(:)));
if nargin<6 % calculating beta_check, theta_check, W
    % fprintf('Calculating W_i..\n');
    W = cell(1,M);
    big_Z = zeros(sum(n), M*q);
    long_Z = zeros(sum(n), q);
    long_X = zeros(sum(n), p);
    long_Y = zeros(sum(n),1);
    G = zeros(sum(n),1);
    for i=1:M
        big_Z(1+sum(n(1:i-1)):sum(n(1:i)), 1+(i-1)*q:i*q) = Z{i};
        long_Z(1+sum(n(1:i-1)):sum(n(1:i)), :) = Z{i};
        long_X(1+sum(n(1:i-1)):sum(n(1:i)), :) = X{i};
        long_Y(1+sum(n(1:i-1)):sum(n(1:i))) = Y{i};
        G(1+sum(n(1:i-1)):sum(n(1:i))) = i;
    end
    lme = fitlmematrix([long_X, big_Z], long_Y, long_Z, G, 'CovariancePattern', 'Isotropic','FitMethod','REML');
    [psi, sigma] = covarianceParameters(lme);
    for i=1:M
        W{i} = (sigma*eye(n(i))+Z{i}*psi{1}*Z{i}')\eye(n(i));
    end
    % fprintf('Initialization done.\n');

    %% Step 1: Calculate check parameters
    % fprintf('Step 1: Calculate check parameters.\n');
    beta_check = zeros(M, p);
    theta_check = zeros(M, q);
    tic;
    for i=1:M
        T = [X{i},Z{i}];
        check = (T'*W{i}*T) \ T'*W{i}*Y{i};
        beta_check(i,:) = check(1:p);
        theta_check(i,:) = check(p+1:end);
    end
    timecost(1) = toc;
    % fprintf('Step 1 done. Timecost: %.6fs\n',timecost(1));
end

%% ADMM
% BIC tuning
min_BIC = Inf;
% fprintf('Step 2: ADMM\n');

% B
% fprintf('Calculating B.\n');
B = zeros(q*M*(M-1)/2, q*M);
t = 0;
for i=1:M
    for j=i+1:M
        t = t+1;
        B(t*q-q+1:t*q, :) = [zeros(q,q*(i-1)), eye(q), zeros(q,q*(j-i-1)), -eye(q), zeros(q,q*(M-j))];
    end
end
B = B';

% big Z
% fprintf('Calculating Z.\n');
bigZ = blkdiag(Z{:});

% big W
% fprintf('Calculating W.\n');
bigW = blkdiag(W{:});

% D
% fprintf('Calculating D.\n');
D = [];
for i=1:M
    D = blkdiag(D, Z{i}'*W{i}*X{i});
end


% bigTheta
Theta_Check = theta_check';
Theta_Check = Theta_Check(:);
bigTheta = Theta_Check;

% beta
Beta_Check = beta_check';
Beta_Check = Beta_Check(:);

% delta, nu
delta = B'*bigTheta;

t=0;
max_lambda = 10*norm(delta)/M/(M-1)*2;
lambda_full = 0:0.2*max_lambda:max_lambda;
lambda_full = lambda_full(2:end);
max_iter = 150;
tic;
for lambda = lambda_full
    t = t+1;
    % fprintf('Lambda: %.4f\n', lambda);

    % bigTheta
    bigTheta = Theta_Check;
    % delta, nu
    delta = B'*bigTheta;
    nu = zeros(q*M*(M-1)/2,1);
    k = 0;
    converge = true;
    while true
        k = k+1;
        %% beta
        LHS = 0;
        RHS = 0;
        for i=1:M
            K_i = X{i}'*W{i}*X{i};
            LHS = LHS + K_i;
            RHS = RHS + K_i * beta_check(i,:)' + X{i}'*W{i}*Z{i}*(theta_check(i,:)'-bigTheta((i-1)*q+1:i*q));
        end
        beta_k = LHS \ RHS;

        %% Theta
        Beta_k = beta_k(:,ones(M,1));
        Beta_k = Beta_k(:);
        bigTheta = (bigZ'*bigW*bigZ + N*(B*B')) \ (bigZ'*bigW*bigZ*Theta_Check + D*(Beta_Check-Beta_k) + N*B*(delta - nu));

        %% delta
        delta = B'*bigTheta + nu;
        for i=1:M*(M-1)/2
            temp = delta((i-1)*q+1:i*q);
    %         if norm(temp) <= 1.85*lambda  % TLP
    %             delta((i-1)*q+1:i*q) = max(1-lambda/norm(temp),0)*temp;
    %         end
            if norm(temp) <= 3.7*lambda % MCP
                delta((i-1)*q+1:i*q) = max(1-lambda/norm(temp),0)*temp/(1-1/3.7);
            end
        end

        %% nu
        nu = nu + B'*bigTheta - delta;

        %% Convergence confirm
        r = B'*bigTheta - delta;
        %fprintf('||r||^2 = %.10f\n', norm(r)^2);
        if norm(r)^2 < epsilon
            break;
        end

        if k > max_iter
    %         fprintf('Doesnt converge.\n');
            fprintf('%d iterations exceeded.\n', max_iter);
            converge = false;
            break;
        end
    end
    if converge
        fprintf('k=%d, converged.\n', k);
        beta_t = beta_k;
        theta_t = reshape(bigTheta, q, M);
        theta_t = theta_t';
        alpha_t = uniquetol(theta_t, 1e-3, 'byrows', true);
        S_t = size(alpha_t,1);
    %     fprintf('Estimated S: %d\n', S_t);
    
        theta_val = zeros(size(theta_check_val));
        subgroup_val = cell(1,S_t);
        for s=1:S_t
            subgroup_val{s} = [];
        end
        [~, theta_val] = estimate_groups(subgroup_val, alpha_t, theta_check_val);
        BIC = bic(X_val, Y_val, Z_val, beta_t', theta_val, S_t);
    %     fprintf('BIC: %.4f\n', BIC);

%         subgroup_t = cell(1,S_t);
%         for s=1:S_t
%             subgroup_t{s} = [];
%         end
%         for i=1:M
%             for s=1:S_t
%                 if sum(abs((theta_t(i,:)-alpha_t(s,:)))>1e-3)==0
%                     subgroup_t{s} = [subgroup_t{s}, i];
%                     break;
%                 end
%             end
%         end
%         subgroup_full{t} = subgroup_t;
    else
        continue;
    end

    if BIC<min_BIC
        min_BIC = BIC;
        beta = beta_t';
        theta = theta_t;
        alpha = alpha_t;
%         best_lambda = lambda;
        S = S_t;
        best_t = t;
    end

end
timecost(2) = toc;
fprintf('Best Lambda: No.%d\n', best_t);
% fprintf('Step 2 done. Time cost: %.6fs\n', timecost(2));



subgroup = cell(1,S);
for s=1:S
    subgroup{s} = [];
end
tic;
for i=1:M
    for s=1:S
        if sum(abs((theta(i,:)-alpha(s,:)))>1e-6)==0
            subgroup{s} = [subgroup{s}, i];
            break;
        end
    end
end
timecost(3) = toc;
timecost = sum(timecost(2:6));
fprintf('Total time cost: %.6fs\n', timecost);
end
