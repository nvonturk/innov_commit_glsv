%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: value_function_iteration_outer.m
% Author: Nicholas von Turkovich
% Date: 10/25/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is a wrapper for the VFI code and performs the additional
% step of ensuring labor market clearing in the event that labor supply is
% inelastic

function [econparams] = value_function_iteration_outer(econparams, x_vector_overload)
    
    %% If labor supply is inelastic
    if ~econparams.elastic_labor_flag
        
        % Iteration counter and distance
        iter = 1;
        w_diff = 10000;

        % Set up bounds on omega
        omega_ceiling = 1;
        omega_floor = 0;
        
        % Check to make sure initial guess within 0,1
        if ~(econparams.omega < 1 || econparams.omega > 0)
            error("Omega must be between 0 and 1 with inelastic labor supply")
        end
        
        % Check to make sure entry is directed
        if econparams.entry_flag == -1
            error("Entry with inelastic labor supply must be directed")
        end
    
        while (abs(w_diff) > econparams.w_tol)
            
            if iter ~= 1
               if (econparams.w > 1)
                   omega_floor = econparams.omega;
                   econparams.omega = (econparams.omega + omega_ceiling)/2;
               elseif (econparams.w < 1)
                   omega_ceiling = econparams.omega;
                   econparams.omega = (econparams.omega + omega_floor)/2;
               end
            end
            
            % Recompute production labor demand
            econparams = compute_labor_demand(econparams);
        
            % Converge given current omega
            econparams = value_function_iteration_inner(econparams);
            
            % Compute the stationary distribution of states
            econparams = compute_stationary_distribution(econparams);

            % Compute the new wage
            econparams = compute_wage(econparams);
            
            % Compute difference from final value
            w_diff = econparams.w - 1;

            % Increment iteration variable
            iter = iter + 1;
        end
    %% If labor supply is elastic
    else
        
        % Check to make sure omega == 1
        if ~(econparams.omega == 1)
            error("Omega must be 1 with elastic labor supply")
        end
        
        % Solve with VFI
        if nargin > 1
            econparams = value_function_iteration_inner(econparams, x_vector_overload);
        else
            econparams = value_function_iteration_inner(econparams);
        end
        
    end

end
