clear;
load('real_data/climate.mat');
load('real_data/unit_GLS.mat');
M = size(X, 2);
n = zeros(1,M);
for i=1:M
    n(i) = size(X{i}, 1);
end
p = size(X{1}, 2);
q = size(Z{1}, 2);

rng(12345679);
% timecost_dishes = zeros(1,10);
% timecost_kmns = zeros(1,10);
% S_dishes = zeros(1,10);
% S_kmns = zeros(1,10);
% error_dishes = zeros(1,10);
% error_kmns = zeros(1,10);
% file_dishes = fopen(sprintf('real_data/Report_DISHES.csv'),'w');
% fprintf(file_dishes, 'Nu,Timecost,S,Error\n');
% file_dishes = fopen(sprintf('real_data/Report_DISHES.csv'),'a');
% 
% file_kmns = fopen(sprintf('real_data/Report_kmeans.csv'), 'w');
% fprintf(file_kmns, 'Timecost,S,Error\n');
% file_kmns = fopen(sprintf('real_data/Report_kmeans.csv'),'a');

for K=2:10

    % train-test split
%     test_size = floor(M*0.2);
%     valid_size = floor(M*0.1);
%     rand_idx = randperm(M);
%     test_idx = rand_idx(1:test_size);
%     valid_idx = rand_idx(test_size+1:test_size+valid_size);
%     train_idx = rand_idx(test_size+valid_size+1:end);

    % DISHES
%     fprintf('DISHES:\n');
    
%     for nu = [0.001, 0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 0.99]
%     fprintf('nu = %.3f:\n', nu);
%     [beta_dishes, alpha_dishes, theta_dishes, subgroup_dishes, timecost_dishes, best_nu] = dishes(X, Z, Y, train_idx, valid_idx, [], theta_U, W, Sigma_big);
%     S_dishes = size(subgroup_dishes,2);
%     subgroup_dishes = recover_full_index(subgroup_dishes, train_idx); % recover the train set index into full data set index
%     [subgroup_dishes, theta_dishes] = estimate_groups(subgroup_dishes, alpha_dishes, theta_U);
%     error_dishes = pred_err(X, Z, Y, beta_dishes, theta_dishes, test_idx);
%     fprintf('S: %d, error: %.4f.\n', S_dishes, error_dishes);
%     fprintf(file_dishes, sprintf('%.3f,%.6f,%.2f,%.4f\n',best_nu,timecost_dishes, S_dishes, error_dishes));
%     end

    % k-means
    fprintf('k-means:\n');
    train_idx = 1:M;
    valid_idx = 1:M;
    [beta_kmns, alpha_kmns, theta_kmns, subgroup, ~, timecost_kmns] = kmeans(X, Z, Y, train_idx, valid_idx, theta_U, W, K);
    S_kmns = size(subgroup,2);
%     subgroup_kmns = recover_full_index(subgroup_kmns, train_idx); % recover the train set index into full data set index
%     [subgroup_kmns, theta_kmns] = estimate_groups(subgroup_kmns, alpha_kmns, theta_U);
    error_kmns = pred_err(X, Z, Y, beta_kmns, theta_kmns);
    fprintf('S: %d, error: %.4f.', K, error_kmns);
    fprintf(file_kmns, sprintf('%.6f,%.2f,%.4f\n',timecost_kmns, S_kmns, error_kmns));
    color = zeros(1,M);
    for s=1:K
        for m=subgroup{s}
            color(m) = s;
        end
    end
    color(48)=[];
    file = fopen('real_color_kmns.txt','a');
    fprintf(file, sprintf('# k=%d\ncolor_exp = [', K));
    for m=1:M-1
        fprintf(file, sprintf('%d,', color(m)));
    end
    fprintf(file, sprintf('%d]\n', color(m)));
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
end