function [x, y, u, er] = bvp2d(a, b, c, d, n, m, k, beta_x_func, beta_y_func, gamma_func, g_func, h_func, q_func, f_func, bound_conditions)

    er = 0;

    # Discretize the domain
    h = [(b - a) / (n - 1), (d - c) / (m - 1)];
    calc_x = @(i) a + (i - 1) * h(1);
    calc_y = @(i) c + (i - 1) * h(2);

    x = arrayfun(calc_x, [1 : n]);
    y = arrayfun(calc_y, [1 : m]);

    x_value = @(i) x(clamp(idivide(int32(i), int32(n), "fix") + 1, 0, n));
    y_value = @(i) y(clamp(mod(i, n) + 1, 0, n));

    # Aproximate derivatives by finite differences
    a_coeff = @(i) gamma_func(x_value(i), y_value(i)) + 2 * k * (1 / (h(1) * h(1)) + 1 / (h(2) * h(2)));
    b_coeff = @(i) (-k / (h(1) * h(1))) - (beta_x_func(x_value(i), y_value(i)) / (2 * h(1)));
    c_coeff = @(i) (-k / (h(1) * h(1))) + (beta_x_func(x_value(i), y_value(i)) / (2 * h(1)));
    d_coeff = @(i) (-k / (h(2) * h(2))) - (beta_y_func(x_value(i), y_value(i)) / (2 * h(2)));
    e_coeff = @(i) (-k / (h(2) * h(2))) + (beta_y_func(x_value(i), y_value(i)) / (2 * h(2)));

    # Create pentadiagonal linear system
    N = n * m;

    A = sparse(N, N);

    A(1,1) = a_coeff(1);
    A(1,2) = c_coeff(1);
    A(1,1+n) = e_coeff(1);

    for I = 2 : n
        A(I,I-1)=b(I);
        A(I,I)=a(I);
        A(I,I+1)=c(I);
        A(I,I+n)=e(I);
    end

    for I = n + 1 : ((m-1)*n)
        A(I,I-n) = d_coeff(I);
        A(I,I-1) = b_coeff(I);
        A(I,I) = a_coeff(I);
        A(I,I+1) = c_coeff(I);
        A(I,I+n) = e_coeff(I);
    end

    for I = (((m-1)*n)+1):((m*n)-1)
        A(I,I-n) = d_coeff(I);
        A(I,I-1) = b_coeff(I);
        A(I,I) = a_coeff(I);
        A(I,I+1) = c_coeff(I);
    end

    A(N,N) = a_coeff(N);
    A(N,N-1) = b_coeff(N);
    A(N,N-1) = d_coeff(N);

    f = zeros(N);
    for i = 1 : n
        for j = 1 : m
            f(i * n + j) = f_func(x(i), y(j));
        end
    end

    # Apply bound conditions
    top_bound = [1 : n];

    left_bound = [1 : m];
    left_bound = arrayfun(@(a) a * n, left_bound);

    bottom_bound = [1 : n];
    bottom_bound = arrayfun(@(a) (m - 1) * n + a, bottom_bound);

    right_bound = [0 : (m - 1)];
    right_bound = arrayfun(@(a) a * n + 1, right_bound);

    for i = 1 : numel(bound_conditions)
        bound_condition = bound_conditions(i);

        point = [x_value(bound_condition.condition_index), y_value(bound_condition.condition_index)];

        if strcmp(bound_condition.condition_type, "value")
            row = zeros(N, 1);
            row(bound_condition.condition_index) = 1;
            A(bound_condition.condition_index,:) = row;
            f(bound_condition.condition_index) = g_func(point(1), point(2));
        elseif strcmp(bound_condition.condition_type, "derivative")
            if find(top_bound == bound_condition.condition_index) != []
                A(i,i) += A(i,i-n);
                f(i) += A(i,i-n) * h(2) * h_func(point(1), point(2)) / k;
                A(i,i-n) = 0;
            elseif find(left_bound == bound_condition.condition_index) != []
                A(i,i) += A(i,i+1);
                f(i) += A(i,i+1) * h(1) * h_func(point(1), point(2)) / k;
                A(i,i+1) = 0;
            elseif find(bottom_bound == bound_condition.condition_index) != []
                A(i,i) += A(i,i+n);
                f(i) += A(i,i+n) * h(2) * h_func(point(1), point(2)) / k;
                A(i,i+n) = 0;
            elseif find(right_bound == bound_condition.condition_index) != []
                A(i,i) += A(i,i-1);
                f(i) += A(i,i-1) * h(1) * h_func(point(1), point(2)) / k;
                A(i,i-1) = 0;
            end
        elseif strcmp(bound_condition.condition_type, "mixed")
            if find(top_bound == bound_condition.condition_index) != []
                A(i,i) += A(i,i-n) * (1 - (h(2) * bound_condition.beta_value) / bound_condition.alpha_value);
                f(i) -= A(i,i-n) * ((h(2) * q_func(point(1), point(2))) / bound_condition.alpha_value) * q_func(point(1), point(2));
                A(i,i-n) = 0;
            elseif find(left_bound == bound_condition.condition_index) != []
                A(i,i) += A(i,i-n) * (1 - (h(1) * bound_condition.beta_value) / bound_condition.alpha_value);
                f(i) -= A(i,i-n) * ((h(1) * q_func(point(1), point(2))) / bound_condition.alpha_value) * q_func(point(1), point(2));
                A(i,i-n) = 0;
            elseif find(bottom_bound == bound_condition.condition_index) != []
                A(i,i) += A(i,i-n) * (1 - (h(2) * bound_condition.beta_value) / bound_condition.alpha_value);
                f(i) -= A(i,i-n) * ((h(2) * q_func(point(1), point(2))) / bound_condition.alpha_value) * q_func(point(1), point(2));
                A(i,i-n) = 0;
            elseif find(right_bound == bound_condition.condition_index) != []
                A(i,i) += A(i,i-n) * (1 - (h(1) * bound_condition.beta_value) / bound_condition.alpha_value);
                f(i) -= A(i,i-n) * ((h(1) * q_func(point(1), point(2))) / bound_condition.alpha_value) * q_func(point(1), point(2));
                A(i,i-n) = 0;
            end
        end
    end

    # Solve linear system
    opts.type = "crout";
    opts.droptol = 1e-4;
    [L, U] = ilu(A, opts);
    u = gmres(A, f, 10, 1e-4, 15, L, U);

endfunction