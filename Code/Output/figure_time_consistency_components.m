%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_time_consistency_components.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_time_consistency_components(econparams, econparams_bgp, consistent, optcons_duple_transition)

    % Eps value (decimal change in variable to represent a small change)
    eps_val = 1e-3;

    % Cross section in time (T - 24 and T - 1)
    t_vals = econparams.interval_lengths(1) + 1 - [24,1]./econparams.frequency;

    % Transitions assuming a slightly higher eta2 and slightly lower eta2
    higher_heta2_tx = compute_transition(econparams, econparams_bgp, [consistent.heta_optcons_duple(1), consistent.heta_optcons_duple(2) * (1 + eps_val)], [econparams.rho, econparams.rho], [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);
    lower_heta2_tx = compute_transition(econparams, econparams_bgp, [consistent.heta_optcons_duple(1), consistent.heta_optcons_duple(2) * (1 - eps_val)], [econparams.rho, econparams.rho], [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);

    % Partial derivative of innovation strategies with respect to a small
    % increase in eta2
    dxs_dheta2 = (higher_heta2_tx.x_vector_path - lower_heta2_tx.x_vector_path)/(2*eps_val*consistent.heta_optcons_duple(2));

    % Change in welfare associated with increase in innovation strategy
    dW_dxs = nan(size(optcons_duple_transition.x_vector_path));

    for t = t_vals

        % Only perform partial derivative for s < sbar since the innovation
        % strategy for the leader in the max gap will be 0 no matter the
        % present/future policy
        parfor s = 1:(2*econparams.n) 

            % Adjustment to increase x_sigma(t)
            tx_temp_high.x_vector_path = optcons_duple_transition.x_vector_path;
            tx_temp_high.x_vector_e_path = optcons_duple_transition.x_vector_e_path;
            tx_temp_high.value_function_path = optcons_duple_transition.value_function_path;

            tx_temp_high.x_vector_path(s,t) = optcons_duple_transition.x_vector_path(s,t) * (1 + eps_val);
            [tx_temp_high.mu_path] = shoot_forward(econparams, econparams_bgp.mu, optcons_duple_transition.heta_path, optcons_duple_transition.rho_path, tx_temp_high.x_vector_path, tx_temp_high.x_vector_e_path, optcons_duple_transition.final_leg_start);
            tx_temp_high =  compute_transition_stats(econparams, tx_temp_high, econparams_bgp.Q0, econparams_bgp.Y0, optcons_duple_transition.heta_path, optcons_duple_transition.rho_path, optcons_duple_transition.final_leg_start);

            % Adjustment to decrease x_sigma(t)
            tx_temp_low.x_vector_path = optcons_duple_transition.x_vector_path;
            tx_temp_low.x_vector_e_path = optcons_duple_transition.x_vector_e_path;
            tx_temp_low.value_function_path = optcons_duple_transition.value_function_path;

            tx_temp_low.x_vector_path(s,t) = optcons_duple_transition.x_vector_path(s,t) * (1 - eps_val);
            [tx_temp_low.mu_path] = shoot_forward(econparams, econparams_bgp.mu, optcons_duple_transition.heta_path, optcons_duple_transition.rho_path, tx_temp_low.x_vector_path, tx_temp_low.x_vector_e_path, optcons_duple_transition.final_leg_start);
            tx_temp_low =  compute_transition_stats(econparams, tx_temp_low, econparams_bgp.Q0, econparams_bgp.Y0, optcons_duple_transition.heta_path, optcons_duple_transition.rho_path, optcons_duple_transition.final_leg_start);

            dW_dxs(s, t) = (tx_temp_high.W - tx_temp_low.W)/(2*eps_val*optcons_duple_transition.x_vector_path(s,t));

        end

    end

    iter = 1;
    patterns = ["-", "--"];

    f = figure("name", "time_consistency_components");

    subplot(1,2,1)

    for t = t_vals

        hold on
        plot(-econparams.n:1:econparams.n, dxs_dheta2(:,t), patterns(iter))
        hold off

        iter = iter + 1;

    end

    xline_prop = xline(0, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xlim([-15,15])
    ylim([-1.75,0])
    xlabel("Firm technology position, $\sigma$")
    title("$\frac{\partial x_{\sigma}(t)}{\partial \eta_2}$")
    legend("t = T - " + string((econparams.interval_lengths(1) + 1 - t_vals)*econparams.frequency), "Location", "SouthWest")
    legend boxoff
    grid off

    iter = 1;

    subplot(1,2,2)

    for t = t_vals

        hold on
        plot(-econparams.n:1:econparams.n, dW_dxs(:,t), patterns(iter))
        hold off

        iter = iter + 1;

    end

    xline_prop = xline(0, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xlim([-15,15])
    xlabel("Firm technology position, $\sigma$")
    title("$\frac{\partial W}{\partial x_{\sigma}(t)}$")
    grid off
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")

end



