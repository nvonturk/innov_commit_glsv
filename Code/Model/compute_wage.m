%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_wage.m
% Author: Nicholas von Turkovich
% Date: 11/8/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = compute_wage(econparams)
    
    % Instantiate wage
    w = 0;
    
    % For each gap, compute contribution to wage
    for s = 0:econparams.n
       
        % See Acemoglu Ackigit (2012) equation 23
        w = w + econparams.mu(s+1)*...
            (econparams.production_labor_demand(econparams.n + 1 + s) + ...
             econparams.production_labor_demand(econparams.n + 1 - s) + ...
            (econparams.x_vector(econparams.n + 1 + s)/econparams.B)^(1/econparams.gamma) + ...
            (econparams.x_vector(econparams.n + 1 - s)/econparams.B)^(1/econparams.gamma) + ...
            (econparams.x_vector_e(econparams.n + 1 - s)/econparams.B_e)^(1/econparams.gamma));
    end
    
    econparams.w = w;
    
end