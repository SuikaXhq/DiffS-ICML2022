function main_full(method, is_demo)
if nargin<2
    is_demo = false;
end
if strcmp(method, 'dishes')
    method = 1;
elseif strcmp(method, 'cd')
    method = 2;
elseif strcmp(method, 'kmeans')
    method = 3;
elseif strcmp(method, 'kmeans-open')
    method = 4;
end
method_number = 4;
method_name = cell(1,method_number);
method_name{1} = 'DISHES';
method_name{2} = 'CD';
method_name{3} = 'kmeans';
method_name{4} = 'kmeans-open';
method_case_idx = cell(1,method_number);
method_case_idx{1} = 1:18;
method_case_idx{2} = [1:3, 7:11, 14:15];
method_case_idx{3} = 1:18;
method_case_idx{4} = 1:18;
if is_demo==true
    for i=1:method_number
        method_case_idx{i} = 0;
    end
end

method_rep = cell(1,method_number);
method_rep{1} = 10;
method_rep{2} = 10;
method_rep{3} = 10;
method_rep{4} = 10;
if is_demo==true
    for i=1:method_number
        method_rep{i} = 4;
    end
end


file = fopen(sprintf('results/Report_%s.csv', method_name{method}),'w');
fprintf(file, 'Case,Timecost,S_mean,S_std,NMI,Perfect_recover,RMSE_beta,RMSE_beta_std,RMSE_theta,RMSE_theta_std,RMSE_alpha,RMSE_alpha_std,Error,Error_std,Replicates\n');
fclose(file);

fprintf('Simulated Data Experiments(%s):\n', method_name{method});
for case_number = method_case_idx{method}
% for case_number = 0  % Demo
S_est_full = zeros(1,method_rep{method});
timecost_full = zeros(1,method_rep{method});
NMI_full = zeros(1,method_rep{method});
perfect_full = zeros(1,method_rep{method});
subgroup_est = cell(1,method_rep{method});
beta_est_full = cell(1,method_rep{method});
alpha_est_full = cell(1,method_rep{method});
theta_est_full = cell(1,method_rep{method});
RMSE_beta_full = zeros(1,method_rep{method});
RMSE_theta_full = zeros(1,method_rep{method});
RMSE_alpha_full = zeros(1,method_rep{method});
error_full = zeros(1,method_rep{method});
S_correct = []; % reps that the method gives correct S

load(sprintf('data/Case%d.mat', case_number));

for j = 1:method_rep{method}
    load(sprintf('data/Case%d_Rep%d_unit_GLS.mat', case_number, j));
	fprintf('Method: %s, Case: %d, Replicate: %d\n', method_name{method}, case_number, j);
%     train_idx_full{j} = 1:50;
    switch method
        case 1
            [beta_est_full{j}, alpha_est_full{j}, theta_est_full{j}, subgroup_est{j}, timecost_full(j), ~] = dishes(X_full{j}, Z_full{j}, Y_full{j}, train_idx_full{j}, valid_idx_full{j}, 0.001, theta_U, W, Sigma_big);
        case 2
            [beta_est_full{j}, alpha_est_full{j}, theta_est_full{j}, subgroup_est{j}, timecost_full(j), ~] = cd_fusion(X_full{j}, Z_full{j}, Y_full{j}, train_idx_full{j}, valid_idx_full{j}, beta_U, theta_U, W);
        case 3
            [beta_est_full{j}, alpha_est_full{j}, theta_est_full{j}, subgroup_est{j}, ~, timecost_full(j)] = kmeans(X_full{j}, Z_full{j}, Y_full{j}, train_idx_full{j}, valid_idx_full{j}, theta_U, W);
        case 4
            [beta_est_full{j}, alpha_est_full{j}, theta_est_full{j}, subgroup_est{j}, ~, timecost_full(j)] = kmeans(X_full{j}, Z_full{j}, Y_full{j}, train_idx_full{j}, valid_idx_full{j}, theta_U, W, size(subgroup_full{j},2));
    end
    subgroup_est{j} = recover_full_index(subgroup_est{j}, train_idx_full{j}); % recover the train set index into full data set index
    [subgroup_est{j}, theta_est_full{j}] = estimate_groups(subgroup_est{j}, alpha_est_full{j}, theta_U);
    S_est_full(j) = size(subgroup_est{j},2);
    [NMI_full(j), perfect_full(j)] = nmi(subgroup_full{j}, subgroup_est{j});
    RMSE_beta_full(j) = rmse(beta_est_full{j}, beta_full{j});
    RMSE_theta_full(j) = rmse(theta_est_full{j}, theta_full{j}, test_idx_full{j});
    if S_est_full(j) == size(alpha_full{j}, 1)
        RMSE_alpha_full(j) = rmse_alpha(alpha_est_full{j}, alpha_full{j});
        S_correct(end+1) = j;
    end
    error_full(j) = pred_err(X_full{j}, Z_full{j}, Y_full{j}, beta_est_full{j}, theta_est_full{j}, test_idx_full{j});
    fprintf('S: %d, NMI: %.4f, error: %.4f.', S_est_full(j), NMI_full(j), error_full(j));
    if perfect_full(j)
        fprintf(' Perfect recovery.');
    end
    fprintf('\n\n');
end

NMI = mean(NMI_full);
timecost = median(timecost_full);
perfect_recover = mean(perfect_full);
S_mean = mean(S_est_full);
S_std = std(S_est_full);
RMSE_beta = mean(RMSE_beta_full);
RMSE_beta_std = std(RMSE_beta_full);
RMSE_theta = mean(RMSE_theta_full);
RMSE_theta_std = std(RMSE_theta_full);
RMSE_alpha_full = RMSE_alpha_full(S_correct);
RMSE_alpha = mean(RMSE_alpha_full);
RMSE_alpha_std = std(RMSE_alpha_full);
error_est = mean(error_full);
error_std = std(error_full);


file = fopen(sprintf('results/Report_%s.csv', method_name{method}),'a');
fprintf(file, sprintf('%d,%.6f,%.2f,%.2f,%.6f,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%d\n', case_number, timecost, S_mean, S_std, NMI, perfect_recover, RMSE_beta, RMSE_beta_std, RMSE_theta, RMSE_theta_std, RMSE_alpha, RMSE_alpha_std, error_est, error_std, method_rep{method}));
fclose(file);
save(sprintf('results/Case%d_%s.mat', case_number, method_name{method}), 'beta_full', 'alpha_full', 'theta_full', 'S_est_full','timecost_full','NMI_full','perfect_full','subgroup_est', '-v7.3');
clear -regexp *_full;
end

end % end function