function BIC = bic(X, Y, Z, beta_est, theta_est, S)
M = size(X,2);
p = size(X{1},2);
q = size(Z{1},2);
n = zeros(M,1);
for i=1:M
    n(i) = size(X{i},1);
end
N = sum(n);

SSE = 0;
for i=1:M
    SSE = SSE + norm(Y{i} - X{i}*beta_est' - Z{i}*theta_est(i,:)')^2;
end

BIC = log(SSE/N) + log(log(p+S*q))*log(N)*(p+S*q)/N;
end