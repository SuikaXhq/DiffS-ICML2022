load('real_data/climate.mat');
% load('real_data/microarray.mat');
M = size(X, 2);
n = zeros(1,M);
for i=1:M
    n(i) = size(X{i}, 1);
end
p = size(X{1}, 2);
q = size(Z{1}, 2);

% train-test split
rng(12345679); % seed
train_idx = cell(1,10);
test_idx = cell(1,10);
valid_idx = cell(1,10);
for i=1:10
    test_size = floor(M*0.2);
    valid_size = floor(M*0.1);
    rand_idx = randperm(M);
    test_idx{i} = rand_idx(1:test_size);
    valid_idx{i} = rand_idx(test_size+1:test_size+valid_size);
    train_idx{i} = rand_idx(test_size+valid_size+1:end);
end

% calculate unit GLS
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
lme = fitlmematrix([long_X, big_Z], long_Y, long_Z, group_var, 'CovariancePattern', 'Isotropic','FitMethod','REML');
[psi, sigma] = covarianceParameters(lme);

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
save(sprintf('real_data/climate_unit_GLS.mat'), 'theta_U', 'beta_U', 'Sigma_big', 'W', 'train_idx', 'valid_idx', 'test_idx', '-v7.3');
% save(sprintf('real_data/microarray_unit_GLS.mat'), 'theta_U', 'beta_U', 'Sigma_big', 'W', 'train_idx', 'valid_idx', 'test_idx', '-v7.3');