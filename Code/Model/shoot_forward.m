%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: shoot_forward.m
% Author: Nicholas von Turkovich
% Date: 11/15/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu_path] = shoot_forward(econparams, mu_start, heta_path, rho_path, x_vector_path, x_vector_e_path, final_leg_start, statedep_cutoff)
    
    %% Setup 
    mu_path = nan(length(mu_start),length(heta_path)+1);
    
    % We place the initial mu distribution in the first column of the
    % transition path; then, the innovation decisions will determine how
    % the distribution evolves over time
    mu_path(:,1) = mu_start;
    
    % Create "temp" econparams object to house the value function and
    % compute innovation rates over the transition
    econparams_temp = econparams;
    
    %% Shoot forward 
    
    % The innovation decisions have been computed and these will be the
    % driving force between changes in mu until we reach a new BGP
    
    % Iterate for length of the transition
        
    for i = 1:length(heta_path)
        
        % Pull out the current eta value and the innovation rates
        if nargin == 8 && ~isempty(statedep_cutoff)
            econparams_temp = update_heta(econparams_temp, heta_path(:,i), statedep_cutoff);
        else
            econparams_temp = update_heta(econparams_temp, heta_path(i));
        end
        econparams_temp.rho = rho_path(i);
        econparams_temp.x_vector = x_vector_path(:,i);
        econparams_temp.x_vector_e = x_vector_e_path(:,i);
        
        % Compute the transition matrix for mu (however if on final leg no
        % need since innovation rates are now constant given the final
        % value function)
        if i <= final_leg_start
            M = compute_mu_transition(econparams_temp);
        end
        
        % Compute the mu for next period
        mu_path(:,i+1) = M * mu_path(:,i);
        
    end

end