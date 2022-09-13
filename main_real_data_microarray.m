clear;
load('real_data/microarray.mat');
load('real_data/microarray_unit_GLS.mat');
X = X_full;
Y = Y_full;
Z = Z_full;
M = size(X, 2);
n = zeros(1,M);
for i=1:M
    n(i) = size(X{i}, 1);
end
p = size(X{1}, 2);
q = size(Z{1}, 2);

rng(12345679);
timecost_dishes = zeros(1,10);
timecost_kmns = zeros(1,10);
timecost_metag = zeros(1,10);
S_dishes = zeros(1,10);
S_kmns = zeros(1,10);
S_metag = zeros(1,10);
error_dishes = zeros(1,10);
error_kmns = zeros(1,10);
error_metag = zeros(1,10);
% file_dishes = fopen(sprintf('real_data/Report_microarray_DISHES.csv'),'w');
% fprintf(file_dishes, 'Nu,Timecost,S,S_std,Error,Error_std\n');
% file_dishes = fopen(sprintf('real_data/Report_microarray_DISHES.csv'),'a');

% file_kmns = fopen(sprintf('real_data/Report_microarray_kmeans.csv'), 'w');
% fprintf(file_kmns, 'Timecost,S,S_std,Error,Error_std\n');
% file_kmns = fopen(sprintf('real_data/Report_microarray_kmeans.csv'),'a');

% file_metag = fopen(sprintf('real_data/Report_microarray_MeTaG.csv'), 'w');
% fprintf(file_kmns, 'Timecost,S,S_std,Error,Error_std\n');
% file_metag = fopen(sprintf('real_data/Report_microarray_MeTaG.csv'),'a');

for i=1:10

    % % train-test split
    % test_size = floor(M*0.2);
    % valid_size = floor(M*0.1);
    % rand_idx = randperm(M);
    % test_idx = rand_idx(1:test_size);
    % valid_idx = rand_idx(test_size+1:test_size+valid_size);
    % train_idx = rand_idx(test_size+valid_size+1:end);

    % DISHES
    fprintf('DISHES:\n');
    
%     for nu = [0.001, 0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 0.99]
%     fprintf('nu = %.3f:\n', nu);
    [beta_dishes, alpha_dishes, theta_dishes, subgroup_dishes, timecost_dishes(i), ~] = dishes(X, Z, Y, train_idx{i}, valid_idx{i}, 0.001, theta_U, W, Sigma_big);
    S_dishes(i) = size(subgroup_dishes,2);
    subgroup_dishes = recover_full_index(subgroup_dishes, train_idx{i}); % recover the train set index into full data set index
    [subgroup_dishes, theta_dishes] = estimate_groups(subgroup_dishes, alpha_dishes, theta_U);
    error_dishes(i) = pred_err(X, Z, Y, beta_dishes, theta_dishes, test_idx{i});
    fprintf('S: %d, error: %.4f.\n', S_dishes(i), error_dishes(i));
    % fprintf(file_dishes, sprintf('%.3f,%.6f,%.2f,%.4f\n',best_nu,timecost_dishes, S_dishes, error_dishes));
%     end

    % k-means
    fprintf('k-means:\n');
    [beta_kmns, alpha_kmns, theta_kmns, subgroup_kmns, ~, timecost_kmns(i)] = kmeans(X, Z, Y, train_idx{i}, valid_idx{i}, theta_U, W);
    S_kmns(i) = size(subgroup_kmns,2);
    subgroup_kmns = recover_full_index(subgroup_kmns, train_idx{i}); % recover the train set index into full data set index
    [subgroup_kmns, theta_kmns] = estimate_groups(subgroup_kmns, alpha_kmns, theta_U);
    error_kmns(i) = pred_err(X, Z, Y, beta_kmns, theta_kmns, test_idx{i});
    fprintf('S: %d, error: %.4f.\n', S_kmns(i), error_kmns(i));
    % fprintf(file_kmns, sprintf('%.6f,%.2f,%.4f\n',timecost_kmns, S_kmns, error_kmns));
    % 
    % 
    % % CD
    % fprintf('CD Fusion:\n');
    % [beta_cd, alpha_cd, theta_cd, subgroup_cd, timecost_cd, ~] = cd_fusion(X, Z, Y, train_idx, valid_idx, beta_U, theta_U, W);
    % S_cd = size(subgroup_cd,2);
    % subgroup_cd = recover_full_index(subgroup_cd, train_idx); % recover the train set index into full data set index
    % [subgroup_cd, theta_cd] = estimate_groups(subgrou_cd, alpha_cd, theta_U);
    % error_cd = pred_err(X, Z, Y, beta_cd, theta_cd, test_idx);
    % fprintf('S: %d, error: %.4f.', S_cd, error_cd);

    % MeTaG
    fprintf('MeTaG:\n');
    [beta_metag, alpha_metag, theta_metag, subgroup_metag, timecost_metag(i)] = MeTaG(X, Z, Y, train_idx{i}, valid_idx{i}, theta_U);
    S_metag(i) = size(subgroup_metag,2);
    subgroup_metag = recover_full_index(subgroup_metag, train_idx{i}); % recover the train set index into full data set index
    [subgroup_metag, theta_metag] = estimate_groups(subgroup_metag, alpha_metag, theta_U);
    error_metag(i) = pred_err(X, Z, Y, beta_metag, theta_metag, test_idx{i});
    fprintf('S: %d, error: %.4f.\n', S_metag(i), error_metag(i));
    % fprintf(file_metag, sprintf('%.6f,%.2f,%.4f\n',timecost_metag, S_metag, error_metag));
end

S_std_dishes = std(S_dishes);
S_std_kmns = std(S_kmns);
S_std_metag = std(S_metag);
S_dishes = mean(S_dishes);
S_kmns = mean(S_kmns);
S_metag = mean(S_metag);
error_std_dishes = std(error_dishes);
error_std_kmns = std(error_kmns);
error_std_metag = std(error_metag);
error_dishes = mean(error_dishes);
error_kmns = mean(error_kmns);
error_metag = mean(error_metag);
timecost_dishes = median(timecost_dishes);
timecost_kmns = median(timecost_kmns);
timecost_metag = median(timecost_metag);

file = fopen(sprintf('real_data/Report_microarray.csv'), 'w');
fprintf(file, 'Timecost,S,S_std,Error,Error_std\n');
file = fopen(sprintf('real_data/Report_microarray.csv'),'a');

fprintf(file, sprintf('%.6f,%.2f,%.2f,%.4f,%.4f\n',timecost_dishes, S_dishes, S_std_dishes, error_dishes, error_std_dishes));
fprintf(file, sprintf('%.6f,%.2f,%.2f,%.4f,%.4f\n',timecost_kmns, S_kmns, S_std_kmns, error_kmns, error_std_kmns));
fprintf(file, sprintf('%.6f,%.2f,%.2f,%.4f,%.4f\n',timecost_metag, S_metag, S_std_metag, error_metag, error_std_metag));