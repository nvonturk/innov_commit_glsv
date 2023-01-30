%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_optimal_transition.m
% Author: Nicholas von Turkovich
% Date: 11/15/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W] = compute_optimal_transition(econparams, econparams_bgp, heta_changes, rho_changes, interval_lengths, statedep_cutoff)
   
    if nargin > 5 && ~isempty(statedep_cutoff)
        [transition] = compute_transition(econparams, econparams_bgp, heta_changes, rho_changes, interval_lengths, [], statedep_cutoff);
    else    
        [transition] = compute_transition(econparams, econparams_bgp, heta_changes, rho_changes, interval_lengths);
    end
    
    % Return the (negative) welfare
    W = -1*transition.W;

end