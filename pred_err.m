function predict_error = pred_err(X, Z, Y, beta, theta, test_idx)
% beta - 1xp vector
% theta - Mxq vector

if nargin < 6
    test_idx = 1:size(X,2);
end
X = X(test_idx);
Z = Z(test_idx);
Y = Y(test_idx);
theta = theta(test_idx, :);

M = size(X,2);
n = zeros(M,1);
for i=1:M
    n(i) = size(X{i},1);
end
N = sum(n);

Y_pred = zeros(N,1);
Y_full = zeros(N,1);
for i=1:M
    offset = sum(n(1:i-1))+1;
    Y_pred(offset:offset+n(i)-1) = X{i}*beta'+Z{i}*theta(i,:)';
    Y_full(offset:offset+n(i)-1) = Y{i};
end


predict_error = mean((Y_full-Y_pred).^2);


end