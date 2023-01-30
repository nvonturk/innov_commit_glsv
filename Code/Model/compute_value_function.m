%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_value_function.m
% Author: Nicholas von Turkovich
% Date: 11/4/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = compute_value_function(econparams, x_vector_overload)
    
    %% Compute innovation decisions
    econparams = compute_innovation_decision(econparams);

    % If the maximum x_vector value is greater than x_bar, adjust x_bar
    if max(econparams.x_vector) > econparams.x_bar
        econparams.x_bar = max(econparams.x_vector)*2;
    end
    
    if max(econparams.x_vector_e) > econparams.x_bar
        econparams.x_bar = max(econparams.x_vector_e)*2;
    end

    %% Compute continuation value

    % Determine probability of each new state given sources of state
    % change: own innovation, competitor innovation, entry, patent expiry
    innovation_probs = (econparams.innovation_transitions .* econparams.x_vector');
    
    if nargin > 1
        innovation_probs_c = (econparams.innovation_transitions .* x_vector_overload');
        
        new_state_probs = innovation_probs + ...
            rot90(innovation_probs_c, 2) + ...
            econparams.patent_transitions .* econparams.heta;
    else
        new_state_probs = innovation_probs + ...
            rot90(innovation_probs, 2) + ...
            econparams.patent_transitions .* econparams.heta;
    end
    
    % Depending on type of entry, add in entry impact on new state probs
    % differently
    if econparams.entry_flag == -1
        new_state_probs = new_state_probs + econparams.entry_impact * econparams.x_vector_e;
    else
        new_state_probs = new_state_probs + econparams.entry_impact .* [econparams.x_vector_e', flipud(econparams.x_vector_e(1:econparams.n))'];
    end

    % Uniformization step where we adjust the probabilities of transition
    new_state_probs = new_state_probs ./ ((2 + abs(econparams.entry_flag))*econparams.x_bar + econparams.heta);
    new_state_probs_e = econparams.x_vector_e ./ ((2 + abs(econparams.entry_flag))*econparams.x_bar + econparams.heta);

    new_state_probs_sum = sum(new_state_probs,1);
    
    % Compute the probability of staying in place and deduct the
    % probability of an entrant unseating 
    for i = 1:size(new_state_probs,2)
        
        new_state_probs(i,i) = 1 - new_state_probs_sum(i);
        
        if i < econparams.n + 1
            if econparams.entry_flag == -1
                new_state_probs(i,i) = new_state_probs(i,i) - new_state_probs_e;
            else
                new_state_probs(i,i) = new_state_probs(i,i) - new_state_probs_e(i);
            end
        elseif i == econparams.n + 1
            if econparams.entry_flag == -1
                new_state_probs(i,i) = new_state_probs(i,i) - 0.5*new_state_probs_e;
            else
                new_state_probs(i,i) = new_state_probs(i,i) - 0.5*new_state_probs_e(i);
            end
        end
        
    end

    % Compute the continuation values
    continuation_values = (econparams.value_function' * new_state_probs)';
    
    %% Compute R&D costs
    rd_costs = ((1 - econparams.tau_rd)*econparams.omega)*((econparams.x_vector./econparams.B).^(1/econparams.gamma));
    
    %% Uniformization step for operating profits and discount factor
    operating_profits = (econparams.profits*(1 - econparams.tau) - rd_costs)/((2 + abs(econparams.entry_flag))*econparams.x_bar + econparams.heta);
    discount_factor = ((2 + abs(econparams.entry_flag))*econparams.x_bar + econparams.heta)/((2 + abs(econparams.entry_flag))*econparams.x_bar + econparams.heta + econparams.rho);
    
    %% Compute the new value function
    econparams.value_function = operating_profits + discount_factor*continuation_values;
end