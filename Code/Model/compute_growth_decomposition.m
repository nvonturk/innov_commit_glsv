%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_growth_decomposition.m
% Author: Nicholas von Turkovich
% Date: 5/16/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [baseline_tx] = compute_growth_decomposition(econparams, econparams_bgp_start, heta_changes, rho_changes, interval_lengths)
    
    % Compute baseline transition
    baseline_tx = compute_transition(econparams, econparams_bgp_start, heta_changes, rho_changes, interval_lengths);
    baseline_tx_high = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2)*1.1], rho_changes, interval_lengths);
    baseline_tx_low = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2)*0.9], rho_changes, interval_lengths);
    baseline_alt_tx_high = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2)*1.1], rho_changes, interval_lengths, baseline_tx.x_vector_path(:,end-1), baseline_tx.x_vector_path, 2);
    baseline_alt_tx_low = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2)*0.9], rho_changes, interval_lengths, baseline_tx.x_vector_path(:,end-1), baseline_tx.x_vector_path, 2);

    figure(1)
    subplot(1,2,1)
    plot(-econparams.n:1:econparams.n, baseline_tx.x_vector_path(:,end-1), 'k-')
    hold on
    plot(-econparams.n:1:econparams.n, baseline_tx_high.x_vector_path(:,end-1), 'r-')
    plot(-econparams.n:1:econparams.n, baseline_alt_tx_high.x_vector_path(:,end-1), 'g-')
    hold off
    xlabel("$\sigma$")
    ylabel("$x_{\sigma}$")
    legend("$\eta_2$", "$\eta_2*1.1$", "$\eta_2*1.1$ (comp. fixed)")
    subplot(1,2,2)
    plot(-econparams.n:1:econparams.n, baseline_tx.x_vector_path(:,end-1), 'k-')
    hold on
    plot(-econparams.n:1:econparams.n, baseline_tx_low.x_vector_path(:,end-1), 'b-')
    plot(-econparams.n:1:econparams.n, baseline_alt_tx_low.x_vector_path(:,end-1), 'g-')
    hold off
    xlabel("$\sigma$")
    ylabel("$x_{\sigma}$")
    legend("$\eta_2$", "$\eta_2*0.9$", "$\eta_2*0.9$ (comp. fixed)")
    
    figure(2)
    subplot(1,2,1)
    plot(-econparams.n:1:econparams.n, baseline_alt_tx_high.x_vector_path(:,end-1) - baseline_alt_tx_high.x_vector_path(:,end-1), 'k-')
    hold on
    plot(-econparams.n:1:econparams.n, baseline_alt_tx_high.x_vector_path(:,end-2) - baseline_alt_tx_high.x_vector_path(:,end-1), 'r-')
    plot(-econparams.n:1:econparams.n, baseline_alt_tx_high.x_vector_path(:,end-3) - baseline_alt_tx_high.x_vector_path(:,end-1), 'g-')
    hold off
    xlabel("$\sigma$")
    ylabel("$x_{\sigma}$")
    legend("$BGP$", "$T-1$", "$T-2$ (comp. fixed)")
    subplot(1,2,2)
    plot(-econparams.n:1:econparams.n, baseline_alt_tx_high.value_function_path(:,end) - baseline_alt_tx_high.value_function_path(:,end) , 'k-')
    hold on
    plot(-econparams.n:1:econparams.n, baseline_alt_tx_high.value_function_path(:,end-1) - baseline_alt_tx_high.value_function_path(:,end) , 'r-')
    plot(-econparams.n:1:econparams.n, baseline_alt_tx_high.value_function_path(:,end-2) - baseline_alt_tx_high.value_function_path(:,end) , 'g-')
    hold off
    xlabel("$\sigma$")
    ylabel("$v_{\sigma}$")
    legend("$BGP$", "$T-1$", "$T-2$ (comp. fixed)")
    
    % Eps value
    eps_val = 1e-11;
    
    higher_tx = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2) + eps_val], rho_changes, interval_lengths);
    lower_tx = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2) - eps_val], rho_changes, interval_lengths);
    
    higher_bgp = value_function_iteration_robust(update_heta(econparams, heta_changes(2) + eps_val), false);
    lower_bgp = value_function_iteration_robust(update_heta(econparams, heta_changes(2) - eps_val), false);
    
    % Calculate derivatives with respect to epsilon
    dg = (higher_tx.g - lower_tx.g)/(2*eps_val);
    dmus = (higher_tx.mu_path - lower_tx.mu_path)/(2*eps_val);
    dxs = (higher_tx.x_vector_path - lower_tx.x_vector_path)/(2*eps_val);
    dxs_elasticity = dxs * eps_val ./ baseline_tx.x_vector_path;
    dxs_bgp = (higher_bgp.x_vector - lower_bgp.x_vector)/(2*eps_val);
    dxs_bgp_elasticity = dxs_bgp * eps_val ./ baseline_tx.x_vector_path;

    % Compute the composition effect
    inc_expected_increase = (econparams.inc_innovation_increase .* baseline_tx.x_vector_path);
    
    composition_effect = log(econparams.lambda)*sum((dmus .* (inc_expected_increase(econparams.n + 1 : (2*econparams.n + 1), :) + inc_expected_increase(econparams.n + 1 : -1: 1, :))), 1);
    
    % Compute the intensive margin
    d_inc_expected_increase = (econparams.inc_innovation_increase .* dxs);
    
    intensive_margin = log(econparams.lambda)*sum((baseline_tx.mu_path .* (d_inc_expected_increase(econparams.n + 1 : (2*econparams.n + 1), :) + d_inc_expected_increase(econparams.n + 1 : -1: 1, :))), 1);
    
    figure(2)
    subplot(1,2,1)
    plot(-econparams.n:1:econparams.n, dxs(:,1/econparams.frequency))
    hold on
    plot(-econparams.n:1:econparams.n, dxs(:,12/econparams.frequency))
    plot(-econparams.n:1:econparams.n, dxs(:,24/econparams.frequency))
    plot(-econparams.n:1:econparams.n, dxs_bgp)
    hold off
    legend(strcat("t = ", num2str(1), " years"), strcat("t = ", num2str(12), " years"), strcat("t = ", num2str(24), " years"), "BGP", "Location", "SouthEast")
    xlabel("$\sigma$")
    ylabel("$\frac{\partial x_{\sigma}(t)}{\partial \epsilon}$")
    
    subplot(1,2,2)
    plot(-econparams.n:1:econparams.n, dxs_elasticity(:,1/econparams.frequency))
    hold on
    plot(-econparams.n:1:econparams.n, dxs_elasticity(:,12/econparams.frequency))
    plot(-econparams.n:1:econparams.n, dxs_elasticity(:,24/econparams.frequency))
    hold off
    legend(strcat("t = ", num2str(1), " years"), strcat("t = ", num2str(12), " years"), strcat("t = ", num2str(24), " years"), "Location", "SouthEast")
    xlabel("$\sigma$")
    ylabel("$\frac{\partial x_{\sigma}(t)}{\partial \epsilon} * \frac{\epsilon}{x_{\sigma}(t)}$")


    
    
    % Compute transition with alternative other firm strategies
        
    tx_vec_high = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2) + eps_val], rho_changes, interval_lengths, baseline_tx.x_vector_path(:,end-1), baseline_tx.x_vector_path, 2);
    tx_vec_low = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2) - eps_val], rho_changes, interval_lengths, baseline_tx.x_vector_path(:,end-1), baseline_tx.x_vector_path, 2);
    
%     dxdsigma = tx_vec_high.
    
    ts = [12, 144, 288];
    
    figure()
    
    for j = 1:length(ts)
        for k = 1:length(sigmas)

            subplot(3,3,3*(j-1) + k)
            hold on
            for i = 1:length(eps_mult)
                plot(eps_mult(i)*eps_val, tx_vec{i}.x_vector_path(sigmas(k),ts(j)), 'ko');
            end
            hold off
            title(strcat("t = ", num2str(ts(j)), ", $\sigma$ = ", num2str(sigmas(k) - econparams.n - 1)))
            xlabel(strcat("$\epsilon$"))
        
        end
    end
    
    higher_tx_kw = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2) + eps_val], rho_changes, interval_lengths, baseline_tx.x_vector_path(:,end-1), baseline_alt_tx.x_vector_path, 2);
    lower_tx_kw = compute_transition(econparams, econparams_bgp_start, [heta_changes(1), heta_changes(2) - eps_val], rho_changes, interval_lengths, baseline_tx.x_vector_path(:,end-1), baseline_alt_tx.x_vector_path, 2);
    dxskw = (higher_tx_kw.x_vector_path - lower_tx_kw.x_vector_path)/(2*eps_val);
    
    d_inc_expected_increase_kw = (econparams.inc_innovation_increase .* dxskw);
    
    kw_effect = log(econparams.lambda)*sum((baseline_alt_tx.mu_path .* (d_inc_expected_increase_kw(econparams.n + 1 : (2*econparams.n + 1), :) + d_inc_expected_increase_kw(econparams.n + 1 : -1: 1, :))), 1);
    
    strategic_effect = log(econparams.lambda)*sum((baseline_tx.mu_path .* (d_inc_expected_increase(econparams.n + 1 : (2*econparams.n + 1), :) + d_inc_expected_increase(econparams.n + 1 : -1: 1, :) -  d_inc_expected_increase_kw(econparams.n + 1 : (2*econparams.n + 1), :) - d_inc_expected_increase_kw(econparams.n + 1 : -1: 1, :) )), 1);
    
    

end