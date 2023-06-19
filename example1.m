function example1()

    a = 0;
    b = 1;
    c = 0;
    d = 1;

    ns = [10, 25, 50, 100];
    ms = [10, 25, 50, 100];

    bound_condition_top.condition_type = "value";
    bound_condition_top.bound = "top";

    bound_condition_right.condition_type = "value";
    bound_condition_right.bound = "right";

    bound_condition_bottom.condition_type = "value";
    bound_condition_bottom.bound = "bottom";

    bound_condition_left.condition_type = "value";
    bound_condition_left.bound = "left";

    bound_conditions = [bound_condition_top, bound_condition_right, bound_condition_bottom, bound_condition_left];

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