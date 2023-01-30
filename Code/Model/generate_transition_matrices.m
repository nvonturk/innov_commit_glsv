%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: generate_transition_matrices.m
% Author: Nicholas von Turkovich
% Date: 12/21/2021
% Note(s): Procuces the transition matrices that dictate transitions in the
% state variable after innovation events or patent expiry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = generate_transition_matrices(econparams)
    
    %% Patent Expiry
    
    % Patent expiry transition matrix
    patent_transition_follower = zeros((econparams.n*2 + 1),econparams.n);
    
    % Zeta controls whether patent expiry completely diffuses knowledge or
    % decreases the gap by a single step
    for column_index = 1:(econparams.n)
        patent_transition_follower(column_index+1,column_index) = (1 - econparams.zeta);
        patent_transition_follower(econparams.n+1,column_index) = patent_transition_follower(econparams.n+1,column_index) + ...
            econparams.zeta;
    end
    
    % Patent transitions have symmetric but opposite effect on leaders
    patent_transition_leader = rot90(patent_transition_follower,2);
    
    % Assemble entire transition matrix; tied industries have no diffusion
    econparams.patent_transitions = [patent_transition_follower, ...
                                     zeros((2*econparams.n+1),1),...
                                     patent_transition_leader];
                                 
    % Create vector of overall expiry rates
    econparams.heta_vector = (sum(econparams.patent_transitions, 1)*econparams.heta)';
    
    %% Incumbent Innovation
    
    % Innovation transition matrix (rows are future states, columns are
    % current states) for incumbents
    econparams.innovation_transitions = zeros((econparams.n*2 + 1));
    
    for column_index = 1:(econparams.n*2 + 1)
       % If follower, with probability phi an innovation leads to a tie/leapfrog;
       % otherwise, advance one step
       if column_index <= (econparams.n)
          econparams.innovation_transitions(column_index+1, column_index) = (1 - econparams.phi);
          econparams.innovation_transitions(econparams.n+1+econparams.leap,column_index) = econparams.innovation_transitions(econparams.n+1+econparams.leap,column_index) + econparams.phi;
       % If tied or leader, advancement is one step
       elseif (column_index > (econparams.n)) && (column_index < (econparams.n*2 + 1))
          econparams.innovation_transitions(column_index+1, column_index) = 1;
       % If maximum gap leader, assume that an innovation leads to leader
       % staying place (hence why innovation rate for this state is 0)
       else
           econparams.innovation_transitions(column_index, column_index) = 1;
       end
    end
    
    %% Entrant Innovation
    
    % Innovation transition matrix for entrants (only of size n+1 since an
    % entrant replaces only the follower/tied firm in each of n+1 gaps)
    econparams.innovation_transitions_e = zeros((econparams.n*2 + 1), (econparams.n + 1));
    
    for column_index = 1:(econparams.n+1)
      % With probability phi_e the entrant
      % experiences drastic innovation, otherwise slow catchup
      if column_index <= econparams.n
          econparams.innovation_transitions_e(column_index+1, column_index) = (1 - econparams.phi_e);
          econparams.innovation_transitions_e(econparams.n + 1 + econparams.leap_e, column_index) = econparams.innovation_transitions_e(econparams.n + 1 + econparams.leap_e, column_index) + econparams.phi_e;
      else
          econparams.innovation_transitions_e(column_index + 1, column_index) = 1;
      end
      
    end
    
    % If no entry then zero out the transition matrix
    econparams.innovation_transitions_e = econparams.innovation_transitions_e * abs(econparams.entry_flag);
    
    % Entry impact matrix for incumbents, this shows incumbents where they
    % end up under entry (note that followers drop out hence why there is
    % no mass for the first n columns; half of tied firms drop one spot)
    econparams.entry_impact = [zeros(2*econparams.n + 1, econparams.n),...
        [zeros(econparams.n-1,1); 0.5; zeros(econparams.n + 1,1)],...
        rot90(econparams.innovation_transitions_e(:,(1:econparams.n)), 2)];
    
    % Growth contribution matrices (upon innovation)
    inc_innovation_increase = zeros(2*econparams.n + 1, 1);
    ent_innovation_increase = zeros(econparams.n + 1, 1);
    
    for s = 1:(econparams.n*2 + 1)
        
        for shat = 1:(econparams.n*2 + 1)
            
            if s <= econparams.n
                num_innovations = max(0,shat - econparams.n - 1);
            else
                num_innovations = max(0,shat - s);
            end
            
            inc_innovation_increase(s) = inc_innovation_increase(s) + num_innovations * econparams.innovation_transitions(shat, s);
            
            if s <= econparams.n + 1
                ent_innovation_increase(s) = ent_innovation_increase(s) + num_innovations * econparams.innovation_transitions_e(shat, s);
            end
            
        end
        
    end    
    
    econparams.inc_innovation_increase = inc_innovation_increase;
    econparams.ent_innovation_increase = ent_innovation_increase;

end