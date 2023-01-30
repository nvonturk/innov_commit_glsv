%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_hjb.m
% Author: Nicholas von Turkovich
% Date: 11/10/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value_function] = compute_hjb(econparams, x_vector_overload, statedep_cutoff)
    
    % Determine probability of each new state given sources of state
    % change: own innovation, competitor innovation, entry, patent expiry
    if nargin > 1 && ~isempty(x_vector_overload)       
        new_state_probs = (econparams.innovation_transitions .* econparams.x_vector') + ...
        rot90(econparams.innovation_transitions .* x_vector_overload', 2) + ...
        econparams.patent_transitions .* econparams.heta;
    elseif nargin > 2 && ~isempty(statedep_cutoff)
        new_state_probs = (econparams.innovation_transitions .* econparams.x_vector') + ...
        rot90(econparams.innovation_transitions .* econparams.x_vector', 2) + ...
        econparams.patent_transitions .* econparams.heta_vector';
    else
        new_state_probs = (econparams.innovation_transitions .* econparams.x_vector') + ...
        rot90(econparams.innovation_transitions .* econparams.x_vector', 2) + ...
        econparams.patent_transitions .* econparams.heta;
    end

    if econparams.entry_flag == -1
        new_state_probs = new_state_probs + econparams.entry_impact * econparams.x_vector_e;
    else
        new_state_probs = new_state_probs + econparams.entry_impact .* [econparams.x_vector_e', flipud(econparams.x_vector_e(1:econparams.n))'];
    end

    if nargin > 1 && ~isempty(x_vector_overload)
        continuation_values = (econparams.value_function' * new_state_probs)' - ...
            (econparams.x_vector .* econparams.value_function) - ...
            (flipud(x_vector_overload) .* econparams.value_function) - ...
            (econparams.heta_vector .* econparams.value_function);
    else
        continuation_values = (econparams.value_function' * new_state_probs)' - ...
            (econparams.x_vector .* econparams.value_function) - ...
            (flipud(econparams.x_vector) .* econparams.value_function) - ...
            (econparams.heta_vector .* econparams.value_function);
    end
    
    if econparams.entry_flag == -1
        continuation_values = continuation_values - econparams.x_vector_e .* econparams.value_function;
    else
        continuation_values = continuation_values - ([econparams.x_vector_e', flipud(econparams.x_vector_e(1:econparams.n))']' .* econparams.value_function);
    end
    
    rd_costs = ((1 - econparams.tau_rd)*econparams.omega)*((econparams.x_vector./econparams.B).^(1/econparams.gamma));
    operating_profits = (econparams.profits*(1 - econparams.tau) - rd_costs);
    
    rhs = operating_profits + continuation_values;

    value_function = rhs/econparams.rho;

end
