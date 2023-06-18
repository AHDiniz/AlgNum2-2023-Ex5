function asymptotic_analysis(a, b, c, d, ns, ms, k, beta_x_func, beta_y_func, gamma_func, g_func, h_func, q_func, f_func, bound_conditions, solution, output_dir)

    results = cell(numel(ns));

    for i = 1 : numel(ns)
    
        n = ns(i);
        m = ms(i);

        [x, y, approx_u, er] = bvp2d(a, b, c, d, n, m, k, beta_x_func, beta_y_func, gamma_func, g_func, h_func, q_func, f_func, bound_conditions);

        real_u = zeros(n * m);
        for j = 1 : n
            for k = 1 : m
                real_u(j * n + k) = solution(x(j), y(k));
            end
        end

        results{i} = {x, y, approx_u, real_u};

        [xx, yy] = meshgrid(x, y);
        
        % hf = figure();
        % surf(xx, yy, approx_u);
        % title(sprintf("Solução Aproximada n = %d, m = %d", n, m));
        % print(hf, sprintf("%s/n_%d-m_%d_approx.png", output_dir, n, m), "-dpng");
        
        % hf = figure();
        % surf(xx, yy, real_u);
        % title(sprintf("Solução Real n = %d, m = %d", n, m));
        % print(hf, sprintf("%s/n_%d-m_%d_real.png", output_dir, n, m), "-dpng");

    end

    space_length = (b - a);
    logh = zeros(numel(ns), 1);
    logE = zeros(numel(ns), 1);
    for i = 1 : numel(ns)
        logh(i) = log(space_length / (ns(i) - 1));
        logE(i) = log(norm(results{i}{3} - results{i}{4}, inf));
    end

    [p] = polyfit(logh, logE, 1);

    px = polyval(p, logh);

    hf = figure();
    plot(logh, logE, "o", logh, px);
    title("Taxa de Convergência 1D");
    print(hf, sprintf("%s/convergence.png", output_dir), "-dpng");

endfunction