%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: shoot_back.m
% Author: Nicholas von Turkovich
% Date: 11/10/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_vector_path, x_vector_e_path, value_function_path] = shoot_back(econparams, value_function_end, heta_path, rho_path, final_leg_start, x_vector_path_overload, statedep_cutoff)
    
    %% Setup 
    x_vector_path = nan(length(econparams.x_vector), length(heta_path)+1);
    x_vector_e_path = nan(length(econparams.x_vector_e), length(heta_path)+1);
    value_function_path = nan(length(value_function_end), length(heta_path)+1);
    
    % Firms are forward looking, in last period of transition they will be
    % looking to the BGP value function to determine their innovation rates
    value_function_path(:,end) = value_function_end;
    
    % Create "temp" econparams object to house the value function and
    % compute innovation rates over the transition
    econparams_temp = econparams;
    
    %% Shoot backward 
    
    % Innovation decisions depend on value function next period; value
    % today is a weighted sum between the continuation values associated
    % with innovation rates and the flow profit today
    
    % Iterate for length of the transition
    for i = size(heta_path,2):-1:1
        
        % Pull out the current eta and rho values
        econparams_temp.value_function = value_function_path(:,i+1);

        if nargin == 7 && ~isempty(statedep_cutoff)
            econparams_temp = update_heta(econparams_temp, heta_path(:,i), statedep_cutoff);
        else
            econparams_temp = update_heta(econparams_temp, heta_path(i));
        end

        econparams_temp.rho = rho_path(i);

        % Compute innovation decisions (note that this involves just the
        % value function and transition matrices + invariant parameters of
        % the model)
        econparams_temp = compute_innovation_decision(econparams_temp);

        % Save the innovation decisions
        x_vector_path(:,i) = econparams_temp.x_vector;
        x_vector_e_path(:,i) = econparams_temp.x_vector_e;

        % Given the current rho and values for eta what is the current
        % period value function
        if nargin == 6 && ~isempty(x_vector_path_overload)
            value_function_path(:,i) = (econparams_temp.rho*compute_hjb(econparams_temp, x_vector_path_overload(:,i)) + ...
            econparams_temp.value_function)./(1 + econparams_temp.rho);
        elseif nargin == 7 && ~isempty(statedep_cutoff)
            value_function_path(:,i) = (econparams_temp.rho*compute_hjb(econparams_temp, [], statedep_cutoff) + ...
            econparams_temp.value_function)./(1 + econparams_temp.rho);
        else
            value_function_path(:,i) = (econparams_temp.rho*compute_hjb(econparams_temp) + ...
            econparams_temp.value_function)./(1 + econparams_temp.rho);
        end

%             % Innovation rates and value function don't change because
%             % there's no policy shift yet
%             x_vector_path(:,i) = x_vector_path(:,i+1);
%             x_vector_e_path(:,i) = x_vector_e_path(:,i+1);
%             value_function_path(:,i) = value_function_path(:,i+1);
            
        
    end

end