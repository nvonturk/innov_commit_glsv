%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_optimal_consistent_duple.m
% Author: Nicholas von Turkovich
% Date: 3/14/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W] = compute_optimal_consistent_duple(econparams, econparams_bgp, heta1, heta2_guess, rho_changes, interval_lengths, heta_lb, heta_ub, options, triple_reopt_filename)
    tic
    % Compute consistent duple
    [heta2_consistent, heta2_max, mu_policy_change, W_schedule] = compute_consistent_duple(econparams, econparams_bgp, heta1, heta2_guess, rho_changes, interval_lengths, heta_lb, heta_ub, options);
    duple_time = toc;
    % Return the (negative) welfare
    W = -1*W_schedule;
    
    % Save consistent duple
    consistent_duple.heta1 = heta1;
    consistent_duple.heta2 = heta2_consistent;
    save("../../Data/Intermediate/" + triple_reopt_filename + "_consistent_duple.mat", "-struct", "consistent_duple");
    fprintf("Duple %2.5f, %2.5f found in %2.2f minutes\n", heta1*100/econparams.frequency, heta2_consistent*100/econparams.frequency, duple_time/60);
end