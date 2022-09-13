function [beta, alpha, theta, subgroup] = GLS(X, Z, Y, train_idx, valid_idx, theta_U, beta_U, W)

X = X(train_idx);
Z = Z(train_idx);
Y = Y(train_idx);
theta_U = theta_U(train_idx, :);
beta_U = beta_U(train_idx, :);
W = W(train_idx);

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

theta = theta_U;
alpha = theta;
subgroup = cell(1,M);
for s=1:M
    subgroup{s} = [s];
end
end