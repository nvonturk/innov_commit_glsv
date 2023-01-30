%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_labor_demand.m
% Author: Nicholas von Turkovich
% Date: 11/8/2021
% Note(s): Compute the demand for production labor in each state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = compute_labor_demand(econparams)
    
    % If kappa is valid, then formula for divying up production labor to
    % followers and leaders
    if econparams.kappa ~= -1
        econparams.production_labor_demand = (1/econparams.omega)*(econparams.profits ./ (econparams.markups - 1));
    else
    % If kappa is a flag, only leaders produce (except in case of tied
    % industries)
        econparams.production_labor_demand = zeros((econparams.n*2 + 1),1);
        econparams.production_labor_demand((econparams.n + 1):(econparams.n*2 + 1)) = (1/econparams.omega)*1./(econparams.markups((econparams.n + 1):(2*econparams.n + 1)));
        econparams.production_labor_demand(econparams.n + 1) = 0.5*econparams.production_labor_demand(econparams.n + 1);
    end
    
end