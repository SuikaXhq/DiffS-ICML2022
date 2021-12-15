function [X, Z, Y, beta_0, alpha_0, theta_0, subgroup] = generate_data(M,S,n,p,q,case_number,rep, sigma_u2, sigma_e2)
if nargin<9
    sigma_u2 = 1; % sigma_u^2
    sigma_e2 = 2; % sigma_epsilon^2
end
if nargin<7
    rep = 100;
end
fprintf('Generating data for case %d, %d replicates.\n', case_number, rep);
X_full = cell(1,rep);
Z_full = cell(1,rep);
Y_full = cell(1,rep);
subgroup_full = cell(1,rep);
beta_full = cell(1,rep);
alpha_full = cell(1,rep);
theta_full = cell(1,rep);
train_idx_full = cell(1,rep);
valid_idx_full = cell(1,rep);
test_idx_full = cell(1,rep);

for k=1:rep
X = cell(1,M);
Z = cell(1,M);
Y = cell(1,M);

%% Effects
beta_0 = rand(1,p)*4 - 2;
beta_full{k} = beta_0;
alpha_01 = (-1:2/(S-1):1)*(S^1.4/2);
if mod(S,2)==1
    alpha_01 = alpha_01 +1;
end
alpha_0 = zeros(S,q);
alpha_0(:,1) = alpha_01;
for j=2:q
    alpha_0(:,j) = alpha_0([end, 1:end-1],j-1);
end
alpha_full{k} = alpha_0;

u = zeros(M, q);
for j=1:M
    u(j,:) = mvnrnd(zeros(1,q), sigma_u2*eye(q));
end

%% Subgroups
theta_0 = zeros(M, q);
M_s = mnrnd(M-S, ones(1,S)/S) + 1;
subgroup = cell(1,S);
offset = 1;
for s=1:S 
    if s==1
        subgroup{s} = 1:M_s(s);
    else
        subgroup{s} = (sum(M_s(1:s-1))+1):(sum(M_s(1:s-1))+M_s(s));
    end
    theta_0(offset:offset+M_s(s)-1, :) = alpha_0(s*ones(M_s(s),1), :);
    offset = offset + M_s(s);
end
theta_full{k} = theta_0;

%% Data
for i=1:M
    %% Design matrix
    Sig_d = 0.3*ones(p+q)+0.7*eye(p+q);
    T = mvnrnd(zeros(1,p+q), Sig_d, n);
    X{i} = T(:, 1:p);
    Z{i} = T(:, p+1:end);
    E = mvnrnd(zeros(1,n), sigma_e2*eye(n));
    
    %% Y_i
    Y{i} = X{i}*beta_0' + Z{i}*(theta_0(i,:)+u(i,:))' + E';
end

X_full{k} = X;
Y_full{k} = Y;
Z_full{k} = Z;
subgroup_full{k} = subgroup;

% train-test split
test_size = floor(M*0.2);
valid_size = floor(M*0.1);
rand_idx = randperm(M);
test_idx_full{k} = rand_idx(1:test_size);
valid_idx_full{k} = rand_idx(test_size+1:test_size+valid_size);
train_idx_full{k} = rand_idx(test_size+valid_size+1:end);


end

%% Save
save(sprintf('data/Case%d.mat', case_number), 'X_full','Z_full','Y_full','subgroup_full', 'beta_full', 'alpha_full', 'theta_full', 'train_idx_full', 'valid_idx_full', 'test_idx_full', 'sigma_u2', 'sigma_e2', '-v7.3');


end
