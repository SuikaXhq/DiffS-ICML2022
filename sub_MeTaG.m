function W = sub_MeTaG(X, Y, W_initial, mu, lambda_h)
    iter_num = 100;
    tau = 1;
    L = 10^4;
    [~, m] = size(Y);
    [~, d] = size(X{1});
    C = zeros([m * (m - 1) / 2 m]);
    C_row_num = 1;

    for i = 1:m

        for j = i + 1:m
            C(C_row_num, i) = 1;
            C(C_row_num, j) = -1;
            C_row_num = C_row_num + 1;
        end

    end

    C = lambda_h * C;
    C = sparse(C);
    A = zeros([m * (m - 1) / 2 d]);
    W = W_initial;

    for iter = 1:iter_num
        W_h = W;
        temp_C_W = C * W_h';

        for i = 1:m * (m - 1) / 2
            A(i, :) = temp_C_W(i, :) / mu;
        end

        for i = 1:m * (m - 1) / 2

            if norm(A(i, :), 2) > 1
                A(i, :) = A(i, :) / norm(A(i, :), 2);
            end

        end

        gradient_f = zeros([d m]);

        for i = 1:m
            gradient_f(:, i) = X{i}' * (X{i} * W_h(:, i) - Y(:, i));
        end

        G_f = gradient_f + A' * C;
        W_h_new = W_h - G_f / L;
        tau_new = 2 / (iter + 3);
        W = W_h_new + (1 - tau) / tau * tau_new * (W_h_new - W_h);
    end
