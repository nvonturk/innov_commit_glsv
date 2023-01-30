%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: strategic_decomposition.m
% Author: Nicholas von Turkovich
% Date: 5/16/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = strategic_decomposition(param_file, transition_type)

    %% Gather data

    % Economy variables
    try
        econparams = load("Data/Intermediate/" + param_file + "_econparams.mat");
    catch
        error(strcat("The econparams variable for '", param_file, "' was not found. Ensure that you have fully run the code under switch baseline prior to running this function."))
    end
    
    % Starting (calibrated) BGP
    econparams_bgp = value_function_iteration_robust(econparams, false);

    % Policy to use
    try
        policy = load("Data/Intermediate/" + param_file + "_" + transition_type + ".mat");
    catch
        error(strcat("The policy variable under '", transition_type, "' for the '", param_file, "' case was not found. Ensure that you have fully run the code under switch baseline prior to running this function."))
    end
        
    if transition_type == "consistent"
        policy = policy.heta_optcons_duple;
    elseif transition_type == "commitment"
        policy = policy.heta_opt_duple;
    else
        error("Specify one of two baseline transitions")
    end

    terminal_bgp = value_function_iteration_robust(update_heta(econparams, policy(length(policy))), false);

    if length(policy) == 3
        interval_lengths = econparams.interval_lengths;
    elseif length(policy) == 2
        interval_lengths = [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))];
    else
        interval_lengths = sum(econparams.interval_lengths);
    end

    transition = compute_transition(econparams, econparams_bgp, policy, repmat(econparams.rho,1,length(policy)), interval_lengths);

    % Eps value
    eps_val = 1e-4;

    %% BGP only

    % Increase/decreate in the patent expiry rate
    terminal_bgp_high = value_function_iteration_robust(update_heta(econparams, (1 + eps_val)*terminal_bgp.heta), false);
    terminal_bgp_low = value_function_iteration_robust(update_heta(econparams, (1 - eps_val)*terminal_bgp.heta), false);

    % Change in growth due to small increase in patent expiry rate
    dg_de = (terminal_bgp_high.g - terminal_bgp_low.g)/(2*eps_val*terminal_bgp.heta);

    % Change in mu due to small increase in patent expiry rate
    dmu_de = (terminal_bgp_high.mu - terminal_bgp_low.mu)/(2*eps_val*terminal_bgp.heta);

    % The "unchanged" innovation strategies at the BGP
    x_vector_overload = terminal_bgp.x_vector;

    % Change in innovation rates due to increase in patent expiry rate holding
    % competitors' strategies constant
    terminal_bgp_high_overload = value_function_iteration_robust(update_heta(econparams, (1 + eps_val)*terminal_bgp.heta), false, x_vector_overload);

    % Change in innovation rates due to decrease in patent expiry rate holding
    % competitors' strategies constant
    terminal_bgp_low_overload = value_function_iteration_robust(update_heta(econparams, (1 - eps_val)*terminal_bgp.heta), false, x_vector_overload);

    % How does this impact the innovation rate of firm in same industry
    dxskw_de = (terminal_bgp_high_overload.x_vector - terminal_bgp_low_overload.x_vector)/(2*eps_val*terminal_bgp.heta);

    % What is the impact when we don't force competitors' strategies to be
    % constant
    dxs_de = (terminal_bgp_high.x_vector - terminal_bgp_low.x_vector)/(2*eps_val*terminal_bgp.heta);

    % The strategic impact on innovation strategies given by the difference
    dxsstrat_de = dxs_de - dxskw_de;

    % Change in innovation rates due to competitor increasing innovation rate
    dxs_dxs = zeros(length(econparams.x_vector), length(econparams.x_vector));

    % Changes in competitors innovation strategies and how they impact other
    % firms
    parfor j = 1:2*econparams.n

        x_vector_overload_high = x_vector_overload;
        x_vector_overload_high(j) = x_vector_overload_high(j)*(1+eps_val);

        x_vector_overload_low = x_vector_overload;
        x_vector_overload_low(j) = x_vector_overload_low(j)*(1-eps_val);

        temp_high = value_function_iteration_robust(update_heta(econparams, terminal_bgp.heta), false, x_vector_overload_high);
        temp_low = value_function_iteration_robust(update_heta(econparams, terminal_bgp.heta), false, x_vector_overload_low);

        dxs_dxs(:,j) = (temp_high.x_vector - temp_low.x_vector)/(2*eps_val*x_vector_overload(j));

    end

    % Determine the contributions due to the escape/trickle effect
    dxs_dxs_escape = zeros(size(dxs_dxs));
    dxs_dxs_trickle = zeros(size(dxs_dxs));

    for i = 1:2*econparams.n + 1

        impacted_idx = i - econparams.n - 1;

        for j = 1:2*econparams.n + 1

            changed_idx = j - econparams.n - 1;

            % Does this contribute to escape competition (i.e., is the
            % impacted_idx a leader and the changed_idx a follower and at the
            % same technology position as the leader or lower gap)
            if impacted_idx >= 0 && changed_idx >= -impacted_idx

                dxs_dxs_escape(i,j) = 1;

            end

            % Does this contribute to trickle down effect (i.e., is the
            % impacted idx a leader and the changed_idx a follower and at a
            % larger gap)
            if impacted_idx >= 0 && changed_idx < -impacted_idx

                dxs_dxs_trickle(i,j) = 1;

            end 

        end

    end

    % Complement gives the other contributions
    dxs_dxs_other = (1 - (dxs_dxs_escape + dxs_dxs_trickle));

    dxs_dxs_escape = dxs_dxs .* dxs_dxs_escape;
    dxs_dxs_trickle = dxs_dxs .* dxs_dxs_trickle;
    dxs_dxs_other = dxs_dxs .* dxs_dxs_other;

    % Contributions of rates/changes in rates of innovation to growth
    inc_expected_increase = (econparams.inc_innovation_increase .* terminal_bgp.x_vector);
    d_inc_expected_increase = (econparams.inc_innovation_increase .* dxs_de); 
    d_inc_expected_increase_kw = (econparams.inc_innovation_increase .* dxskw_de); 
    d_inc_expected_increase_strat = (econparams.inc_innovation_increase .* (dxs_de - dxskw_de));
    d_inc_expected_increase_strat2 = (econparams.inc_innovation_increase .* (dxs_dxs * dxs_de));
    d_inc_expected_increase_strat_escape = (econparams.inc_innovation_increase .* (dxs_dxs_escape * dxs_de));
    d_inc_expected_increase_strat_trickle = (econparams.inc_innovation_increase .* (dxs_dxs_trickle * dxs_de));
    d_inc_expected_increase_strat_other = (econparams.inc_innovation_increase .* (dxs_dxs_other * dxs_de));

    % Breakdown of growth into two broader effects (extensive and intensive
    % margin)
    composition_effect = (dmu_de' * (inc_expected_increase(econparams.n + 1 : -1 : 1) + inc_expected_increase(econparams.n + 1 : 2*econparams.n + 1)));
    intensive_margin = (terminal_bgp.mu' * (d_inc_expected_increase(econparams.n +  1 : -1 : 1) + d_inc_expected_increase(econparams.n + 1 : 2*econparams.n + 1)));

    % Breakdown of intensive margin into the Kremer-Williams effect / strategic
    % channel
    kw_effect = (terminal_bgp.mu' * (d_inc_expected_increase_kw(econparams.n + 1 : -1 : 1) + d_inc_expected_increase_kw(econparams.n + 1 : 2*econparams.n + 1)));
    strategic_effect = (terminal_bgp.mu' * (d_inc_expected_increase_strat2(econparams.n + 1 : -1 : 1) + d_inc_expected_increase_strat2(econparams.n + 1 : 2*econparams.n + 1)));

    % Breakdown of strategic channel into three components (note the other has
    % no impact on growth)
    escape_effect = (terminal_bgp.mu' * (d_inc_expected_increase_strat_escape(econparams.n + 1 : -1 : 1) + d_inc_expected_increase_strat_escape(econparams.n + 1 : 2*econparams.n + 1)));
    trickle_effect = (terminal_bgp.mu' * (d_inc_expected_increase_strat_trickle(econparams.n + 1 : -1 : 1) + d_inc_expected_increase_strat_trickle(econparams.n + 1 : 2*econparams.n + 1)));
    other_effect = (terminal_bgp.mu' * (d_inc_expected_increase_strat_other(econparams.n + 1 : -1 : 1) + d_inc_expected_increase_strat_other(econparams.n + 1 : 2*econparams.n + 1)));

    % Checks
    (dg_de - log(econparams.lambda)*(composition_effect + intensive_margin))/dg_de
    (kw_effect + strategic_effect - intensive_margin)/intensive_margin
    (dg_de - log(econparams.lambda)*(composition_effect + kw_effect + strategic_effect))/dg_de
    (dg_de - log(econparams.lambda)*(composition_effect + kw_effect + escape_effect + trickle_effect + other_effect))/dg_de

    %% With transition

    % Transition varying the policy rate
    policy_high = policy;
    policy_high(end) = policy_high(end)*(1+eps_val);
    policy_low = policy;
    policy_low(end) = policy_low(end)*(1-eps_val);

    tx_high = compute_transition(econparams, econparams_bgp, policy_high, repmat(econparams.rho,1,length(policy)), interval_lengths);
    tx_low = compute_transition(econparams, econparams_bgp, policy_low, repmat(econparams.rho,1,length(policy)), interval_lengths);

    % Change in growth
    dg_de_t = (tx_high.g - tx_low.g)/(2*eps_val*terminal_bgp.heta);

    % Change in mu
    dmu_de_t = (tx_high.mu_path - tx_low.mu_path)/(2*eps_val*terminal_bgp.heta);

    % Change in strategy associated with increase in terminal patent expiry
    % rate
    dxs_de_t = (tx_high.x_vector_path - tx_low.x_vector_path)/(2*eps_val*terminal_bgp.heta);

    % Additional column for the BGP response
    dxs_de_t(:,end) = dxs_de_t(:,end-1);

    % Strategies with no change (baseline)
    x_vector_path_overload = transition.x_vector_path;

    % Create an additional column for solving the BGP change
    x_vector_path_overload(:,end) = x_vector_path_overload(:,end-1);

    % Transition varying the policy rate but holding competitor strategies
    % fixed
    tx_high_overload = compute_transition(econparams, econparams_bgp, policy_high, repmat(econparams.rho,1,length(policy)), interval_lengths, x_vector_path_overload);
    tx_low_overload = compute_transition(econparams, econparams_bgp, policy_low, repmat(econparams.rho,1,length(policy)), interval_lengths, x_vector_path_overload);

    % Change in strategy associated with increase in patent epxiry rate holding
    % competitor strategies constant
    dxskw_de_t = (tx_high_overload.x_vector_path - tx_low_overload.x_vector_path)/(2*eps_val*terminal_bgp.heta);

    % Change in innovation rates due to competitor increasing innovation rate
    % (note for efficiency purposes these will be weighted by the dxs_de_t
    % while aggregating)
    dxs_dxs_t = zeros(size(x_vector_path_overload));
    dxs_dxs_t_unweighted = cell(2*econparams.n + 1, 1);
    for i = 1:2*econparams.n + 1
        dxs_dxs_t_unweighted{i} = dxs_dxs_t;
    end
    dxs_dxs_t_escape = dxs_dxs_t;
    dxs_dxs_t_trickle = dxs_dxs_t;
    dxs_dxs_t_other = dxs_dxs_t;

    for t = length(x_vector_path_overload):-1:1
        
        parfor j = 1:2*econparams.n

            % Small increase in the innovation rate of the competitor in
            % position j at time t
            x_vector_path_overload_high = x_vector_path_overload;
            x_vector_path_overload_high(j,t) = x_vector_path_overload_high(j,t)*(1+eps_val);

            % Symmetric decrease
            x_vector_path_overload_low = x_vector_path_overload;
            x_vector_path_overload_low(j,t) = x_vector_path_overload_low(j,t)*(1-eps_val);

            % Transitions under point changes in strategies at (j,t) without
            % changing the terminal expiry rate
            temp_high = compute_transition(econparams, econparams_bgp, policy, repmat(econparams.rho,1,length(policy)), interval_lengths, x_vector_path_overload_high);
            temp_low = compute_transition(econparams, econparams_bgp, policy, repmat(econparams.rho,1,length(policy)), interval_lengths, x_vector_path_overload_low);

            % Impact on strategies associated with change at (j,t)
            impact = (temp_high.x_vector_path - temp_low.x_vector_path)/(2*eps_val*x_vector_path_overload(j,t))*dxs_de_t(j,t);
            dxs_dxs_t = dxs_dxs_t + impact;
            dxs_dxs_t_unweighted{j} = dxs_dxs_t_unweighted{j} + (temp_high.x_vector_path - temp_low.x_vector_path)/(2*eps_val*x_vector_path_overload(j,t));

            % Portion that impacts the strategies of laggards
            impact_other = impact;
            impact_other(econparams.n+1:2*econparams.n+1,:) =  impact_other(econparams.n+1:2*econparams.n+1,:)*0;
            dxs_dxs_t_other = dxs_dxs_t_other + impact_other;

            % Portion that impacts the strategies of tied/leaders
            impact_growth = impact;
            impact_growth(1:econparams.n,:) = impact_growth(1:econparams.n,:)*0;

            % Determine the indices of the affected firms which, given the
            % value for j, contributed to the trickle and escape effects
            affected_sigma_pr = j - econparams.n - 1;
            affected_trickle = 0:(-1*affected_sigma_pr - 1);
            affected_escape = max(-1*affected_sigma_pr,0):1:econparams.n;

            % Compute the trickle effect component
            impact_growth_trickle = impact_growth;
            impact_growth_trickle(affected_escape + econparams.n + 1,:) =  impact_growth_trickle(affected_escape + econparams.n + 1,:)*0;
            dxs_dxs_t_trickle = dxs_dxs_t_trickle + impact_growth_trickle;

            % Compute the escape effect component
            impact_growth_escape = impact_growth;
            impact_growth_escape(affected_trickle + econparams.n + 1,:) =  impact_growth_escape(affected_trickle + econparams.n + 1,:)*0;
            dxs_dxs_t_escape = dxs_dxs_t_escape + impact_growth_escape;

        end
        fprintf("%2.2f%% done with strategic_decomposition.m\n", 100*(1 + length(x_vector_path_overload) - t)/length(x_vector_path_overload));

    end

    unweighted = zeros(size(dxs_dxs_t));
    for i = 1:length(dxs_dxs_t_unweighted)
        unweighted = unweighted + dxs_dxs_t_unweighted{i};
    end

    % Multiply by the impact that the affected firms have on growth
    inc_expected_increase_t = (econparams.inc_innovation_increase .* transition.x_vector_path);
    d_inc_expected_increase_t = (econparams.inc_innovation_increase .* dxs_de_t); 
    d_inc_expected_increase_kw_t = (econparams.inc_innovation_increase .* dxskw_de_t); 
    d_inc_expected_increase_strat_t = (econparams.inc_innovation_increase .* (dxs_de_t - dxskw_de_t));
    d_inc_expected_increase_strat2_t = (econparams.inc_innovation_increase .* (dxs_dxs_t));
    d_inc_expected_increase_strat_escape_t = (econparams.inc_innovation_increase .* (dxs_dxs_t_escape));
    d_inc_expected_increase_strat_trickle_t = (econparams.inc_innovation_increase .* (dxs_dxs_t_trickle));
    d_inc_expected_increase_strat_other_t = (econparams.inc_innovation_increase .* (dxs_dxs_t_other));

    % Compute the time series for the various effects
    composition_effect_t = sum((dmu_de_t .* (inc_expected_increase_t(econparams.n + 1 : -1 : 1, :) + inc_expected_increase_t(econparams.n + 1 : 2*econparams.n + 1, :))),1);
    intensive_margin_t = sum((transition.mu_path .* (d_inc_expected_increase_t(econparams.n + 1 : -1 : 1, :) + d_inc_expected_increase_t(econparams.n + 1 : 2*econparams.n + 1, :))), 1);
    kw_effect_t = sum((transition.mu_path .* (d_inc_expected_increase_kw_t(econparams.n + 1 : -1 : 1, :) + d_inc_expected_increase_kw_t(econparams.n + 1 : 2*econparams.n + 1, :))), 1);
    strategic_effect_t = sum((transition.mu_path .* (d_inc_expected_increase_strat_t(econparams.n + 1 : -1 : 1, :) + d_inc_expected_increase_strat_t(econparams.n + 1 : 2*econparams.n + 1, :))), 1);
    strategic_effect_escape_t = sum((transition.mu_path .* (d_inc_expected_increase_strat_escape_t(econparams.n + 1 : -1 : 1, :) + d_inc_expected_increase_strat_escape_t(econparams.n + 1 : 2*econparams.n + 1, :))), 1);
    strategic_effect_trickle_t = sum((transition.mu_path .* (d_inc_expected_increase_strat_trickle_t(econparams.n + 1 : -1 : 1, :) + d_inc_expected_increase_strat_trickle_t(econparams.n + 1 : 2*econparams.n + 1, :))), 1);
    strategic_effect_other_t = sum((transition.mu_path .* (d_inc_expected_increase_strat_other_t(econparams.n + 1 : -1 : 1, :) + d_inc_expected_increase_strat_other_t(econparams.n + 1 : 2*econparams.n + 1, :))), 1);

    % Checks
    max(abs(log(econparams.lambda)*(composition_effect_t + intensive_margin_t) - dg_de_t)./abs(dg_de_t))
    max(abs(kw_effect_t + strategic_effect_t - intensive_margin_t)./abs(intensive_margin_t))
    max(abs(strategic_effect_t - strategic_effect_escape_t - strategic_effect_trickle_t)./abs(strategic_effect_t))

    % Save data structures for plotting
    save("Data/Intermediate/" + param_file + "_" + transition_type + "_stratdecomp" + ".mat")

end

