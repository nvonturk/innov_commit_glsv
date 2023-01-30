%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_stationary_distribution.m
% Author: Nicholas von Turkovich
% Date: 10/25/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = compute_stationary_distribution(econparams)
    
    %% Extract data for quicker operations
    n = econparams.n;
    x_vector = econparams.x_vector;
    x_vector_e = econparams.x_vector_e;
    innovation_transitions = econparams.innovation_transitions;
    innovation_transitions_e = econparams.innovation_transitions_e;
    heta_vector = econparams.heta_vector;
    patent_transitions = econparams.patent_transitions;
    entry_flag = econparams.entry_flag;

    %% Equation 1 (density should sum to 1)
    
    % Sum across the density must equal 1
    eq_density = ones(1, n + 1);
    
    %% Equation 2 (outflow and inflows for sigma = 0)
    
    % Outflow of tied industries should equal inflow from other states
    eq_0 = zeros(1, n + 1);
    
    % Outflow due to two incumbents innovating and the possibility of entry
    % but not from patent expiry
    eq_0(1) = -(2*x_vector(n + 1));
    
    if entry_flag == -1
        eq_0(1) = eq_0(1) - x_vector_e;
    else
        eq_0(1) = eq_0(1) - x_vector_e(n+1);
    end
    
    % Inflow coming in from all states apart from 0
    for s_in = 1:n
        
       % Flow out of state s into state 0 due to ...
       eq_0(s_in + 1) = x_vector(n + 1 + s_in)*innovation_transitions(n + 1, n + 1 + s_in) + ... % leader innovation
           x_vector(n + 1 - s_in)*innovation_transitions(n + 1, n + 1 - s_in) + ... % follower innovation
           heta_vector(n + 1 + s_in)*patent_transitions(n + 1, n + 1 + s_in); % patent expiry
       
       if entry_flag == -1
           eq_0(s_in + 1) = eq_0(s_in + 1) + x_vector_e*innovation_transitions_e(n + 1, n + 1 - s_in); % undirected entry
       else
           eq_0(s_in + 1) = eq_0(s_in + 1) + x_vector_e(n + 1 - s_in)*innovation_transitions_e(n + 1, n + 1 - s_in); % entrant innovation
       end
       
    end
    
    %% Equations 3 (inflow and outflow for sigma = 1 to n-1)
    
    % For outflow from states 1 to (n-1) should equal inflows from all
    % other states
    eq_other = zeros((n - 1), (n + 1));
    
    % Outer loop pins down the outflow state; inner loop adds in all the
    % inflows to the outflow state
    for s_out = 1:(n - 1)
        
        % Outflow from state s_out given by ...
        eq_other(s_out, s_out + 1) = -(x_vector(n + 1 + s_out) + ... % leaders innovating out of that gap
            x_vector(n + 1 - s_out)*(1 - innovation_transitions(n + 1 + s_out, n + 1 - s_out)) + ... % followers innovating but not to leader position
            heta_vector(n + 1 + s_out)); % patent expiry
        
        if entry_flag == -1
            eq_other(s_out, s_out + 1) = eq_other(s_out, s_out + 1) - x_vector_e*(1 - innovation_transitions_e(n + 1 + s_out, n + 1 - s_out)); % undirected
        else
            eq_other(s_out, s_out + 1) = eq_other(s_out, s_out + 1) - x_vector_e(n + 1 - s_out)*(1 - innovation_transitions_e(n + 1 + s_out, n + 1 - s_out)); % entrants innovating but not to leader position
        end
        
        % Inflow into state s_out given by ...
        for s_in = 0:(n)
            
            % If on outflow state, skip
            if s_in == s_out
                continue
            end
            
            eq_other(s_out, s_in + 1) = x_vector(n + 1 + s_in)*innovation_transitions(n + 1 + s_out, n + 1 + s_in) + ... % leaders innovating into gap
                x_vector(n + 1 - s_in)*(innovation_transitions(n + 1 + s_out, n + 1 - s_in) + innovation_transitions(n + 1 - s_out, n + 1 - s_in)) + ... % followers innovating into gap as follower/leader
                heta_vector(n + 1 + s_in)*(patent_transitions(n + 1 + s_out, n + 1 + s_in)); % patent expiry sends industry to gap
            
            if entry_flag == -1
                eq_other(s_out, s_in + 1) = eq_other(s_out, s_in + 1) + x_vector_e*(innovation_transitions_e(n + 1 + s_out, n + 1 - s_in) + innovation_transitions_e(n + 1 - s_out, n + 1 - s_in)); % undirected entry
            else
                eq_other(s_out, s_in + 1) = eq_other(s_out, s_in + 1) + x_vector_e(n + 1 - s_in)*(innovation_transitions_e(n + 1 + s_out, n + 1 - s_in) + innovation_transitions_e(n + 1 - s_out, n + 1 - s_in)); % entrant innovating into gap as follower/leader
            end
            
        end
    end
    
    %% Solve for stationary distribution with n+1 equations
    
    % Create A matrix and then solve for mu
    A = [eq_density; eq_0; eq_other];
    b = zeros(n + 1,1);
    b(1) = 1;
    econparams.mu = A \ b;

end