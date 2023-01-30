%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_innovation_decision.m
% Author: Nicholas von Turkovich
% Date: 11/4/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = compute_innovation_decision(econparams)
    
    %% For incumbents

    % Compute the change in value upon innovation
    capital_gain = (econparams.value_function' * econparams.innovation_transitions) -...
        econparams.value_function';
    
    % Determine vector of innovation intensities for incumbents
    % given potential capital cain
    x_vector = econparams.B*((capital_gain' * econparams.gamma * econparams.B) / (econparams.omega * (1 - econparams.tau_rd))).^(econparams.gamma/(1 - econparams.gamma));
        
    % Zero out entries if imaginary
    x_vector(imag(x_vector) ~= 0) = 0;
    
    % No negative rates
    econparams.x_vector = max(x_vector, 0);
    

    
    %% For entrants
    
    if econparams.entry_flag == -1
        
        % Compute the change in value upon innovation
        capital_gain = (econparams.value_function' * econparams.innovation_transitions_e);
        
        % Weight the capital gains by the distribution mu
        capital_gain = sum(capital_gain' .* flipud(econparams.mu));
        
        % Determine vector of innovation intensities for entrants
        % given potential capital cain
        x_vector_e = econparams.B_e*((capital_gain' * econparams.gamma * econparams.B_e) / (econparams.omega * (1 - econparams.tau_rd))).^(econparams.gamma/(1 - econparams.gamma));
        
        % Zero out entries if imaginary
        x_vector_e(imag(x_vector_e) ~= 0) = 0;

        % Zero out vector if not entry
        x_vector_e = abs(econparams.entry_flag) * x_vector_e;

        % No negative rates
        econparams.x_vector_e = max(x_vector_e, 0);
        
        
    else
        % Compute the change in value upon innovation
        capital_gain = (econparams.value_function' * econparams.innovation_transitions_e);

        % Determine vector of innovation intensities for entrants
        % given potential capital cain
        x_vector_e = econparams.B_e*((capital_gain' * econparams.gamma * econparams.B_e) / (econparams.omega * (1 - econparams.tau_rd))).^(econparams.gamma/(1 - econparams.gamma));

        % Zero out entries if imaginary
        x_vector_e(imag(x_vector_e) ~= 0) = 0;

        % Zero out vector if not entry
        x_vector_e = econparams.entry_flag * x_vector_e;

        % No negative rates
        econparams.x_vector_e = max(x_vector_e, 0);
    end
end