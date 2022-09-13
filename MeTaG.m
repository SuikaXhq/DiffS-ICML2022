function [beta_best, alpha_best, theta_best, subgroup_best, timecost] = MeTaG(X, Z, Y, train_idx, valid_idx, theta_U)
    tic;
    H = 2;
    M = size(X,2);
    M_val = size(valid_idx,2);
    p = size(X{1},2);
    q = size(Z{1},2);
    n = size(X{1},1);
    X_full = cell(1,M);
    theta_val = theta_U(valid_idx, :);
    for m=1:M
        X_full{m} = [X{m}, Z{m}];
    end
    X_Train = X_full(train_idx);
    X_Validation = X_full(valid_idx);

    Y_Train = zeros(n, size(train_idx,2));
    Y_Validation = zeros(n, size(valid_idx,2));
    for i=1:size(train_idx,2)
        Y_Train(:, i) = Y{train_idx(i)};
    end
    for i=1:size(valid_idx,2)
        Y_Validation(:, i) = Y{valid_idx(i)};
    end

    Lambda1_Vec = [10,20,50,100,1000];

    [n, T] = size(Y_Train);
    MSE = zeros([1 length(Lambda1_Vec)]);
    W = cell([1 length(Lambda1_Vec)]);
    alpha = cell([1 length(Lambda1_Vec)]);
    beta = cell([1 length(Lambda1_Vec)]);
    theta = cell([1 length(Lambda1_Vec)]);
    subgroup = cell([1 length(Lambda1_Vec)]);
    W_Hierarchy = cell([1 length(Lambda1_Vec)]);

    for i = 1:length(Lambda1_Vec)
        [W{i}, W_Hierarchy{i}] = MeTaG_fun1(X_Train, Y_Train, Lambda1_Vec(i), H);
        beta{i} = W{i}(1:p, :)';
        beta{i} = mean(beta{i}, 1);
        theta{i} = W{i}(p+1:end, :)';
        % theta{i}
        alpha{i} = uniquetol(theta{i}, 1e-3, 'byrows', true);
        S_temp = size(alpha{i}, 1);
        subgroup_temp = cell(1,S_temp);
        for s=1:S_temp
            subgroup_temp{s} = [];
        end
        [subgroup{i}, ~] = estimate_groups(subgroup_temp, alpha{i}, theta{i});
        % Validate
        subgroup_val = cell(1,S_temp);
        for s=1:S_temp
            subgroup_val{s} = [];
        end
        [~, theta_val_temp] = estimate_groups(subgroup_val, alpha{i}, theta_val);
        W_val = [beta{i}(ones(1,M_val), :), theta_val_temp]';
        temp_mse = zeros([1 M_val]);
        for t = 1:M_val
            temp_mse(t) = norm(Y_Validation(:, t) - X_Validation{t} * W_val(:, t), 2)^2 / n;
        end

        MSE(i) = mean(temp_mse);
    end

    [~, I] = min(MSE);
    lambda1 = Lambda1_Vec(I);
    beta_best = beta{I};
    alpha_best = alpha{I};
    theta_best = theta{I};
    subgroup_best = subgroup{I};
    Best_W = W{I};
    Best_W_H = W_Hierarchy{I};
    timecost = toc;
    return;
