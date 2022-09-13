function [W, W_Hierarchy] = MeTaG_fun1(X, Y, Lamb, H)
    [~, m] = size(Y);
    [~, d] = size(X{1});
    mu = 0.01;
    lambda = zeros([1 H]);
    lambda(1) = Lamb;

    for h = 2:H
        lambda(h) = lambda(h - 1) / 1.2;
    end

    W_Hierarchy = cell([1 H]);
    W_0 = zeros([d m]);

    for h = 1:H

        if h == 1
            W_initial = W_0;
            Y_new = Y;
        else
            W_exist = zeros([d m]);

            for k = 1:h - 1
                W_exist = W_exist + W_Hierarchy{k};
            end

            for i = 1:m
                Y_new(:, i) = Y(:, i) - X{i} * W_exist(:, i);
            end

            W_initial = zeros([d m]);

            for i = 1:m
                W_initial(:, i) = lasso(X{i}, Y_new(:, i), 'Lambda', 0.1);
            end

        end

        W_Hierarchy{h} = sub_MeTaG(X, Y_new, W_initial, mu, lambda(h));
    end

    W = zeros([d m]);

    for h = 1:H
        W = W + W_Hierarchy{h};
    end
