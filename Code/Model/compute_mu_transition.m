%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_mu_transition.m
% Author: Nicholas von Turkovich
% Date: 11/15/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = compute_mu_transition(econparams)
    
    %% Setup
    
    % Take out data from econparams for faster computation
    innovation_transitions = econparams.innovation_transitions;
    innovation_transitions_e = econparams.innovation_transitions_e;
    x_vector = econparams.x_vector;
    x_vector_e = econparams.x_vector_e;
    heta_vector = econparams.heta_vector;
    patent_transitions = econparams.patent_transitions;
    n = econparams.n;
    
    % Square matrix to send mu to new distribution
    M = zeros(n + 1);
    
    for shat = 0:n
                
        % Each row gives the weights of the current mu distribution in
        % contributing to the corresponding entry in the resulting
        % distribution
        
        leader_index = n + 1 + shat;
        follower_index = n + 1 - shat;
        
        % First we compute the outflow based leader innovation in tied
        % position
        M(shat+1,shat+1) = (1 - innovation_transitions(leader_index, leader_index) - innovation_transitions(follower_index, leader_index))*x_vector(leader_index);

        % Compute the outflow based on follower innovation
        M(shat+1,shat+1) = M(shat+1,shat+1) + (1 - innovation_transitions(leader_index, follower_index) - innovation_transitions(follower_index, follower_index))*x_vector(follower_index);
        
        % Compute outflow based on entry
        M(shat+1,shat+1) = M(shat+1,shat+1) + (1 - innovation_transitions_e(leader_index, follower_index) - innovation_transitions_e(follower_index, follower_index))*x_vector_e(follower_index);
        
        % Compute outflow due to patent expiry
        M(shat+1,shat+1) = M(shat+1,shat+1) + heta_vector(leader_index);
        
        % Negate the outflow by convention
        M(shat+1,shat+1) = -1*M(shat+1,shat+1);
        
        
        % Then for each other mass in mu, we compute its contribution to
        % inflow to position s
        
        for s = 0:(length(M)-1)
           
            if s == shat
                continue
            end
            
            leader_inflow_index = n + 1 + s;
            follower_inflow_index = n + 1 - s;
            
            % Contributions due to leader innovation
            M(shat+1,s+1) = (innovation_transitions(leader_index, leader_inflow_index) + (shat ~= 0)*innovation_transitions(follower_index, leader_inflow_index))*x_vector(leader_inflow_index);
            
            % Contributions due to follower innovation
            M(shat+1,s+1) = M(shat+1,s+1) + (innovation_transitions(leader_index, follower_inflow_index) + (shat ~= 0)*innovation_transitions(follower_index, follower_inflow_index))*x_vector(follower_inflow_index);

            % Contributions due to entry
            M(shat+1,s+1) = M(shat+1,s+1) + (innovation_transitions_e(leader_index, follower_inflow_index) + (shat ~= 0)*innovation_transitions_e(follower_index, follower_inflow_index))*x_vector_e(follower_inflow_index);
            
            % Contributions due to patent expiry
            M(shat+1,s+1) = M(shat+1,s+1) + (patent_transitions(leader_index, leader_inflow_index) + (shat ~= 0)*patent_transitions(follower_index, follower_inflow_index))*heta_vector(leader_inflow_index);
            
        end
        
        
    end
    
    % Check to make sure that the contributions of each state sum to 1
    if max(abs(sum(M,1))) > 1e-10
        error("Transition matrix is incorrect")
    end
    
    % Add in the identity matrix for immediate multiplication
    M = eye(n + 1) + M;
    

end