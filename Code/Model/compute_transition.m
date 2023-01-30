%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_transition.m
% Author: Nicholas von Turkovich
% Date: 11/10/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [transition] = compute_transition(econparams, econparams_bgp_start, heta_changes, rho_changes, interval_lengths, x_vector_path_overload, statedep_cutoff)
    
    %% Setup path and final BGP
    if econparams.omega ~= 1
        error("Transition only applies for elastic labor supply case")
    end

    % Check to make sure same length for heta_changes and rho_changes
    if size(heta_changes, 2) ~= length(rho_changes)
        error("Vectors for transition must be the same length")
    end

    % Create paths for parameters
    heta_path = nan(size(heta_changes, 1), sum(interval_lengths));
    rho_path = nan(1, sum(interval_lengths));
    
    % For loop to fill in parameters
    start_index = 1;
    
    for i = 1:length(interval_lengths)
        end_index = start_index + (interval_lengths(i)-1);
        heta_path(:, start_index:end_index) = repmat(heta_changes(:, i), 1, length(start_index:end_index));
        rho_path(start_index:end_index) = rho_changes(i);
        start_index = end_index + 1;
    end
    
    % After this loop, find the last start_
    final_leg_start = start_index - interval_lengths(end);
    
    % At time \infty, we are in the BGP associated with the final set of
    % parameters in heta_changes and rho_changes
    econparams_bgp_end = update_heta(econparams, heta_changes(1, end));
    econparams_bgp_end.rho = rho_changes(end);
    
    %% Shoot back
  
    % Starting from the final BGP, we step backwards to determine the
    % innovation decisions for the firms
    
    if nargin == 6 && ~isempty(x_vector_path_overload)
        
        econparams_bgp_end = value_function_iteration_robust(econparams_bgp_end, false, x_vector_path_overload(:,end));
        
        [transition.x_vector_path, transition.x_vector_e_path, transition.value_function_path] = shoot_back(econparams, ...
            econparams_bgp_end.value_function, heta_path, rho_path, final_leg_start, x_vector_path_overload);
        
    elseif nargin == 7 && ~isempty(statedep_cutoff)
        
        econparams_bgp_end = value_function_iteration_robust(econparams_bgp_end, false);
        
        [transition.x_vector_path, transition.x_vector_e_path, transition.value_function_path] = shoot_back(econparams, ...
            econparams_bgp_end.value_function, heta_path, rho_path, final_leg_start, [], statedep_cutoff);
        
        
    else
    
        econparams_bgp_end = value_function_iteration_robust(econparams_bgp_end, false);

        [transition.x_vector_path, transition.x_vector_e_path, transition.value_function_path] = shoot_back(econparams, ...
            econparams_bgp_end.value_function, heta_path, rho_path, final_leg_start);
        
    end
    
    %% Shoot forward
    % Starting from the initial mu
    if nargin == 7 && ~isempty(statedep_cutoff)
        [transition.mu_path] = shoot_forward(econparams, econparams_bgp_start.mu, heta_path, rho_path, transition.x_vector_path, transition.x_vector_e_path, final_leg_start, statedep_cutoff);
    else
        [transition.mu_path] = shoot_forward(econparams, econparams_bgp_start.mu, heta_path, rho_path, transition.x_vector_path, transition.x_vector_e_path, final_leg_start);
    end

    %% Statistics Along the Transition Path
        
    if nargin == 7 && ~isempty(statedep_cutoff)
        transition = compute_transition_stats(econparams, transition, econparams_bgp_start.Q0, econparams_bgp_start.Y0, heta_path, rho_path, final_leg_start, statedep_cutoff);
    else
        transition = compute_transition_stats(econparams, transition, econparams_bgp_start.Q0, econparams_bgp_start.Y0, heta_path, rho_path, final_leg_start);
    end
    
    if abs(transition.W_prod + transition.W_muplab - transition.W) > 1e-5
        error("Check welfare decomposition")
    end
    
    if nargin == 7 && ~isempty(statedep_cutoff)
        check_event = transition.x_vector_path(econparams.n + 1 : 2*econparams.n + 1, 1:end-1) + transition.x_vector_path(econparams.n + 1:-1:1, 1:end-1) + ...
        transition.x_vector_e_path(econparams.n+1:-1:1,1:end-1) + [repmat(transition.heta_path(2,:), length(1:econparams.n + 1 - 1 - statedep_cutoff), 1); repmat(transition.heta_path(1,:), statedep_cutoff, 1); repmat(transition.heta_path(2,:)*0, 1, 1)];
    else
        check_event = transition.x_vector_path(econparams.n + 1 : 2*econparams.n + 1, 1:end-1) + transition.x_vector_path(econparams.n + 1:-1:1, 1:end-1) + ...
        transition.x_vector_e_path(econparams.n+1:-1:1,1:end-1) + transition.heta_path;
    end
    
    if (max(check_event, [], 'all') >= 1)
        error("Frequent events")
    end
    
    if econparams.disp_flag
        disp("Transition completed")
    end

end