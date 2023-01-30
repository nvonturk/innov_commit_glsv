%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: find_consistent_heta.m
% Author: Nicholas von Turkovich
% Date: 11/17/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist, heta2_max, mu_policy_change, W_schedule] = find_consistent_heta(econparams, econparams_bgp, heta1, heta2_guess, rho_changes, interval_lengths, heta_lb, heta_ub, options, statedep_cutoff)
    
    % Given our guess, solve for the transition
    if nargin == 10 && ~isempty(statedep_cutoff)
        transition = compute_transition(econparams, econparams_bgp, [heta1, [heta2_guess; heta2_guess]], rho_changes, interval_lengths, [], statedep_cutoff);
    else
        transition = compute_transition(econparams, econparams_bgp, [heta1, heta2_guess], rho_changes, interval_lengths);
    end
        
    mu_policy_change = transition.mu_path(:,interval_lengths(1) + 1);
    
    % Create the starting economy for the second social planner's problem
    econparams_start = econparams;
    econparams_start.mu = mu_policy_change;
    econparams_start.Q0 = transition.Q(interval_lengths(1) + 1);
    econparams_start.Y0 = transition.C(interval_lengths(1) + 1);
    
    % Solve for the optimal static eta2 conditional on the state of the economy
    % (mu)
    optimal_transition_anon = @(heta2) compute_optimal_transition(econparams, econparams_start, heta2, rho_changes(2), interval_lengths(2));
    [heta2_max, ~] = fmincon(@(heta2) optimal_transition_anon(heta2), [heta2_guess], [], [], [], [], [heta_lb], [heta_ub], [], options);
    
    % Dist > 0 if our guess was low
    dist = heta2_max - heta2_guess;
        
    % Welfare of transition from initial policy maker's perspective
    W_schedule = transition.W;
    
end
