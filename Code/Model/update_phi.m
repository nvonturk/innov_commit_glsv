%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: update_phi.m
% Author: Nicholas von Turkovich
% Date: 11/18/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = update_phi(econparams, phi)

    econparams.phi = phi;
    
    % Innovation transition matrix (rows are future states, columns are
    % current states) for incumbents
    econparams.innovation_transitions = zeros((econparams.n*2 + 1));
    
    for column_index = 1:(econparams.n*2 + 1)
       
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

end