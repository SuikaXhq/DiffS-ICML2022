function [beta, alpha, theta, subgroup, timecost] = CMTL(X, Z, Y, train_idx, valid_idx, theta_U, beta_U, W) 
addpath(genpath('../jiayuzhou-MALSAR-fb97515/MALSAR'));
tic;
if nargin >= 5
    X_val = X(valid_idx);
    Z_val = Z(valid_idx);
    Y_val = Y(valid_idx);
    X = X(train_idx);
    Z = Z(train_idx);
    Y = Y(train_idx);
    if nargin >= 7
        beta_val = beta_U(valid_idx, :);
        theta_U_val = theta_U(valid_idx, :);
        W_val = W(valid_idx);

        beta_U = beta_U(train_idx, :);
        theta_U = theta_U(train_idx,:);
        W = W(train_idx);
    end
end

M = size(X,2);
p = size(X{1},2);
q = size(Z{1},2);
n = zeros(M,1);

lhs = zeros(p,p);
rhs = zeros(p,1);
for i=1:M
    temp = X{i}' * W{i} * X{i};
    lhs = lhs + temp;
    rhs = rhs + temp * beta_U(i, :)';
end
beta = lhs \ rhs;
beta = beta';

X_CMTL = Z;
Y_CMTL = cell(size(Y));
for i=1:M
    Y_CMTL{i} = Y{i} - X{i}*beta';
end

opts.init = 0;
opts.tFlag = 1;
opts.tol = 10^-6;
opts.maxIter = 1000;
rho1 = 10;
rho2 = 10^-1;
% k = 3;
min_BIC = +Inf;
for k=2:10
    fprintf('Test on k = %d.\n', k);
    W_learn = Least_CMTL(X_CMTL, Y_CMTL, rho1, rho2, k, opts);
    W_learn = W_learn';
    alpha_k = uniquetol(W_learn, 1e-3, 'byrows', true);
    S_k = size(alpha_k,1);
    
    % W_learn

    theta_val = zeros(size(theta_U_val));
    subgroup_val = cell(1,S_k);
    for s=1:S_k
        subgroup_val{s} = [];
    end
    [~, theta_val] = estimate_groups(subgroup_val, alpha_k, theta_U_val);
    BIC = bic(X_val, Y_val, Z_val, beta, theta_val, S_k);
    fprintf('Final Estimate for S: %d, BIC: %.4f.\n', S_k, BIC);
    if BIC < min_BIC
        alpha = alpha_k;
        theta = W_learn;
        min_BIC = BIC;
        S = S_k;
        best_k = k;
        [subgroup, ~] = estimate_groups(subgroup_val, alpha_k, W_learn);
    end
end
fprintf('Best estimate: BIC=%.4f, k=%d.\n', min_BIC, best_k)
timecost = toc;

end