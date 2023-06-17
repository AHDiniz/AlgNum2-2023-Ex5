function example1()

    a = 0;
    b = 1;
    c = 0;
    d = 1;

    ns = [10, 25, 50, 100];
    ms = [10, 25, 50, 100];

    bound_conditions = [];

    for i = 1 : numel(ns)
        n = ns(i);
        m = ms(i);

        bound_conditions_array = [];

        top_bound = [1 : n];
        for j = 1 : numel(top_bound)
            bound_condition.condition_type = "value";
            bound_condition.condition_index = top_bound(j);
            bound_conditions_array = [bound_conditions_array, bound_condition];
        end

        left_bound = [1 : m];
        left_bound = arrayfun(@(a) a * n, left_bound);
        for j = 1 : numel(left_bound)
            bound_condition.condition_type = "value";
            bound_condition.condition_index = left_bound(j);
            bound_conditions_array = [bound_conditions_array, bound_condition];
        end

        bottom_bound = [1 : n];
        bottom_bound = arrayfun(@(a) (m - 1) * n + a, bottom_bound);
        for j = 1 : numel(bottom_bound)
            bound_condition.condition_type = "value";
            bound_condition.condition_index = bottom_bound(j);
            bound_conditions_array = [bound_conditions_array, bound_condition];
        end

        right_bound = [0 : (m - 1)];
        right_bound = arrayfun(@(a) a * n + 1, right_bound);
        for j = 1 : numel(right_bound)
            bound_condition.condition_type = "value";
            bound_condition.condition_index = right_bound(j);
            bound_conditions_array = [bound_conditions_array, bound_condition];
        end

        bound_conditions = [bound_conditions, bound_conditions_array];
    end

    k = 1;

    beta_x_func = @(x, y) 1;
    beta_y_func = @(x, y) 20 * y;
    gamma_func = @(x, y) 1;

    g_func = @(x, y) 0;
    h_func = @(x, y) 0;
    q_func = @(x, y) 0;

    f_func = @(x, y) 2.5 * exp(x ^ (9/2)) * (4 * x * x * (41 * y * y - 21 * y - 2) + (18 * (x - 1) * x ^ (9/2) - 9 * (-9 * x ^ (9/2) + 9 * x ^ (11/2) + 15 * x - 11) * x ^(7/2) - 12) * (y - 1) * y - 4 * x * (39 * y * y - 19 * y - 2));

    solution = @(x, y) 10 * x * y * (1 - x) * (1 - y) * exp(x * x * x * x * sqrt(x));

    output_dir = "out/example1";
    
    asymptotic_analysis(a, b, c, d, ns, ms, k, beta_x_func, beta_y_func, gamma_func, g_func, h_func, q_func, f_func, bound_conditions, solution, output_dir);

endfunction