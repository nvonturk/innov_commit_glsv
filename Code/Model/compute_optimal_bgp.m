%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_optimal_bgp.m
% Author: Nicholas von Turkovich
% Date: 11/8/2021
% Note(s): Simply updates the econparams, solves for the BGP, and returns
% the BGP welfare (inverted) for use by fmincon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W] = compute_optimal_bgp(econparams, heta)
    
    % Add in new eta values to economy
    econparams = update_heta(econparams, heta);

    % Find the BGP
    [econparams] = value_function_iteration_robust(econparams, false);
    
    % Return the (negative) welfare
    W = -1*econparams.W;

end