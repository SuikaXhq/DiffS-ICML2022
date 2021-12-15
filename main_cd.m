clear;
file = fopen('results/Report_cd.csv','w');
fprintf(file, 'Case,Timecost,S_mean,S_std,NMI,Perfect_recover,RMSE_beta,RMSE_beta_std,RMSE_theta,RMSE_theta_std,Error,Error_std,Replicates\n');
fclose(file);

fprintf('Simulated Data Experiments(CD):\n');
for case_number = [1:3, 7:15]
% for case_number = 0  % Demo
S_est_full = zeros(1,20);
timecost_full = zeros(1,20);
NMI_full = zeros(1,20);
perfect_full = zeros(1,20);
subgroup_est = cell(1,20);
beta_est_full = cell(1,20);
alpha_est_full = cell(1,20);
theta_est_full = cell(1,20);
RMSE_beta_full = zeros(1,20);
RMSE_theta_full = zeros(1,20);
error_full = zeros(1,20);

load(sprintf('data/Case%d.mat', case_number));
for j = 1:20
    load(sprintf('data/Case%d_Rep%d_unit_GLS.mat', case_number, j));
    fprintf('Method: CD, Case: %d, Replicate: %d\n', case_number, j);
    [beta_est_full{j}, alpha_est_full{j}, theta_est_full{j}, subgroup_est{j}, timecost_full(j), ~] = cd_fusion(X_full{j}, Z_full{j}, Y_full{j}, beta_U, theta_U, W);
    S_est_full(j) = size(subgroup_est{j},2);
    [NMI_full(j), perfect_full(j)] = nmi(subgroup_full{j}, subgroup_est{j});
    RMSE_beta_full(j) = rmse(beta_est_full{j}, beta_full{j});
    RMSE_theta_full(j) = rmse(theta_est_full{j}, theta_full{j});
    error_full(j) = pred_err(X_full{j}, Z_full{j}, Y_full{j}, beta_est_full{j}, theta_est_full{j});
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
error_est = mean(error_full);
error_std = std(error_full);

file = fopen('results/Report_cd.csv','a');
fprintf(file, sprintf('%d,%.6f,%.2f,%.2f,%.6f,%.4f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%d\n', case_number, timecost, S_mean, S_std, NMI, perfect_recover, RMSE_beta, RMSE_beta_std, RMSE_theta, RMSE_theta_std, error_est, error_std, 20));
fclose(file);
save(sprintf('results/Case%d_cd.mat', case_number), 'RMSE_beta_full', 'RMSE_theta_full', 'error_full', 'beta_est_full', 'alpha_est_full', 'theta_est_full', 'S_est_full','timecost_full','NMI_full','perfect_full','subgroup_est', '-v7.3');
clear -regexp *_full;
end
