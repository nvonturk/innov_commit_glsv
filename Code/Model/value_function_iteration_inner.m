%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: value_function_iteration_inner.m
% Author: Nicholas von Turkovich
% Date: 10/25/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function uses value function iteration with uniformization to return
% BGP value functions and investment success rates consistent with 
% discrete-time functional equations for the value function, using the method
% of Acemoglu and Akcigit (2012)'s Online Appendix E.4. 

function [econparams] = value_function_iteration_inner(econparams, x_vector_overload)
    
    %% Initialize loop
    
    iter = 1;
    vdiff = 10000;
    
    % Store old value function
    value_function_old = econparams.value_function;
    
    % For speedy calculations, pull out the objects that do not change on
    % each iteration
    
    %% Loop
    
    while (abs(vdiff) > econparams.v_tol)
        
        % Compute new innovation decisions and value function
        if nargin > 1
            [econparams] = compute_value_function(econparams, x_vector_overload);
        else 
            [econparams] = compute_value_function(econparams);
        end
        
        % If we have undirected entry we will need to update mu as it
        % enters into the decision of the potential entrants
        if econparams.entry_flag == -1
            econparams = compute_stationary_distribution(econparams);
        end

        % Compute difference from prior iteration
        vdiff = max(abs(econparams.value_function - value_function_old));
        
        % Store new value_function as old
        value_function_old = econparams.value_function;
        
        % Increment iteration variable
        iter = iter + 1;
    end

end
