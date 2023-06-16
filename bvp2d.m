function [x, y, u, er] = bvp2d(a, b, c, d, n, m, k, beta_x_func, beta_y_func, gamma_func, g_func, h_func, f_func, bound_conditions)

    er = 0;

    # Discretize the domain
    h = [(b - a) / (n - 1), (d - c) / (m - 1)];
    calc_x = @(i) a + (i - 1) * h(0);
    calc_y = @(i) c + (i - 1) * h(1);

    x = arrayfun(calc_x, [1 : n]);
    y = arrayfun(calc_y, [1 : m]);

    # Aproximate derivatives by finite differences
    a_coeff = @(i) gamma_func(i) + 2 * k * (1 / (h(0) * h(0)) + 1 / (h(1) * h(1)));
    b_coeff = @(i) (-k / (h(0) * h(0))) - (beta_x_func(i) / (2 * h(0)));
    c_coeff = @(i) (-k / (h(0) * h(0))) + (beta_x_func(i) / (2 * h(0)));
    d_coeff = @(i) (-k / (h(1) * h(1))) - (beta_y_func(i) / (2 * h(1)));
    e_coeff = @(i) (-k / (h(1) * h(1))) + (beta_y_func(i) / (2 * h(1)));

    # Create pentadiagonal linear system
    N = n * m;

    A = sparse(N, N);
    for i = 1 : N
        A(i,i) = a_coeff(i);        
    end
    for i = 2 : N
        A(i,i-1) = b_coeff(i);
    end
    for i = 1 : N - 1
        A(i,i+1) = c_coeff(i);
    end
    for i = n + 1 : N
        A(i,i-n) = d_coeff(i);
    end
    for i = 1 : (m - 1) * n
        A(i,i+n) = e_coeff(i);
    end

    f = arrayfun(f_func, [1 : N]);

    # Apply bound conditions
    for i = 1 : numel(bound_conditions)
        bound_condition = bound_conditions(i);

        switch bound_condition.condition_type
            case "value"
                row = zeros(N, 1);
                row(bound_condition.condition_index) = 1;
                A(bound_condition.condition_index) = row;
                f(bound_condition.condition_index) = g_func(bound_condition.condition_index);
            case "derivative"
                indices_a = [1 : n];

                indices_b = [1 : m];
                indices_b = arrayfun(@(a) a * n, indices_b);

                indices_c = [1 : n];
                indices_c = arrayfun(@(a) (m - 1) * n + a, indices_c);

                indices_d = [0 : (m - 1)];
                indices_d = arrayfun(@(a) a * n + 1, indices_d);

                if bound_condition.condition_index == indices_a
                    A(i,i) += A(i,i-n);
                    f(i) += A(i,i-n) * h(1) * h_func(bound_condition.condition_index) / k;
                    A(i,i-n) = 0;
                else if bound_condition.condition_index == indices_b
                    A(i,i) += A(i,i+1);
                    f(i) += A(i,i+1) * h(0) * h_func(bound_condition.condition_index) / k;
                    A(i,i+1) = 0;
                else if bound_condition.condition_index == indices_c
                    A(i,i) += A(i,i+n);
                    f(i) += A(i,i+n) * h(1) * h_func(bound_condition.condition_index) / k;
                    A(i,i+n) = 0;
                else if bound_condition.condition_index == indices_d
                    A(i,i) += A(i,i-1);
                    f(i) += A(i,i-1) * h(0) * h_func(bound_condition.condition_index) / k;
                    A(i,i-1) = 0;
                end
            case "mixed"
            otherwise
        end
    end

    # Solve linear system

endfunction