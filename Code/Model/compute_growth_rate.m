%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_growth_rate.m
% Author: Nicholas von Turkovich
% Date: 12/22/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = compute_growth_rate(econparams)
    
    % Calculate the expected increases to innovation
    inc_expected_increase = econparams.inc_innovation_increase .* econparams.x_vector;
    ent_expected_increase = econparams.ent_innovation_increase .* econparams.x_vector_e;
    
    % Weight increases by density in gaps
    econparams.g = log(econparams.lambda)*(sum(econparams.mu .* (inc_expected_increase(econparams.n + 1:(2*econparams.n + 1)) + inc_expected_increase(econparams.n + 1:-1:1) + ent_expected_increase(econparams.n + 1:-1:1))));
    
    
end