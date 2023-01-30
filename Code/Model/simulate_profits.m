%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: simulate_profits.m
% Author: Nicholas von Turkovich
% Date: 5/4/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = simulate_profits(econparams, transition, sigma_start, t, draws)

    % Starting at the beginning of period t, what is the path that a firm
    % that innovated in t-1 versus one where in t-1 no event happened
    
    % If entry flag is 1, throw error
    if econparams.entry_flag ~= 0
        error("Simulation not set up for entry")
    end
    
    if length(unique(econparams.patent_transitions)) ~= 2
        error("Adjust for varied patent expiry")
    end
    
    % Sum of events
    sum_events = transition.x_vector_path(econparams.n + 1:(2*econparams.n+1), 1:end-1) + ...
        transition.x_vector_path(econparams.n + 1:-1:1, 1:end-1) + ...
        transition.heta_path;
    
    if max(sum_events, [], 'all') > 1
        error("Transition path has period in which event guaranteed to occur")
    end
    
    % Subset the draws
    draws = draws(:, t:end);
    
    % Number of simulations (individual histories)
    N = size(draws, 1);
    
    % Number of time periods
    T = size(draws, 2);
        
    % Instantiate data structures
    firm_state = nan(N, length(sigma_start), T);
    firm_state(:,:,1) = repmat(sigma_start, N, 1);    
    firm_profits = nan(size(firm_state));
    
    % Function to compute R&D costs and profits
    rd_costs = @(x) ((1 - econparams.tau_rd)*econparams.omega)*((x./econparams.B).^(1/econparams.gamma));
    operating_profits = @(sigma, rd) (econparams.profits(econparams.n + 1 + sigma)*(1 - econparams.tau) - rd);
    
    % For each step in time
    for i = 1:T-1
        
        firm_state_curr = squeeze(firm_state(:,:,i));
        firm_state_next = squeeze(firm_state(:,:,i+1));
    
        % For loop over each simulation
        parfor n = 1:N
            
            % Time indicates where in the transition we are (offset by t)
            time = i + (t - 1);

            % Pull out the innovation rates at this point in the transition
            x_vector_temp = transition.x_vector_path(:,time);

            % Pull out the patent expiry rate at this point in the
            % transition
            heta_temp = transition.heta_path(time);
            
            % Compute the strategies of the firms in the simulation based
            % on their current sigma
            firm_strategies = x_vector_temp(econparams.n + 1 + firm_state_curr(n, :));
            competitor_strategies = x_vector_temp(econparams.n + 1 - firm_state_curr(n, :));
            
            % Compute profits
            firm_profits(n,:,i) = operating_profits(firm_state_curr(n, :), rd_costs(firm_strategies));
            
            draw = draws(n,i);
            
            % Construct the intervals
            intervals = [zeros(length(firm_strategies),1), cumsum([firm_strategies, competitor_strategies, repmat(heta_temp, length(firm_strategies), 1)], 2), ones(length(firm_strategies),1)];
            
            % Find index of end of interval for draw
            [~, interval_end] = max(draw - intervals < 0, [], 2);
                
            idx_interval_end = sub2ind(size(intervals), (1:size(intervals,1))', interval_end);
            idx_interval_start = sub2ind(size(intervals), (1:size(intervals,1))', interval_end-1);

            % Calculate what share of the interval it is in
            share_of_interval = (draw - intervals(idx_interval_start))./(intervals(idx_interval_end) - intervals(idx_interval_start));
            
            firm_innovation_destinations = cumsum(econparams.innovation_transitions(:, econparams.n + 1 + firm_state_curr(n, :)), 1);
            [~, firm_innovation_destinations] = max(share_of_interval' - firm_innovation_destinations < 0, [], 1);
            firm_innovation_destinations = firm_innovation_destinations - econparams.n - 1;
            
            competitor_innovation_destinations = cumsum(econparams.innovation_transitions(:, econparams.n + 1 - firm_state_curr(n, :)), 1);
            [~, competitor_innovation_destinations] = max(share_of_interval' - competitor_innovation_destinations < 0, [], 1);
            competitor_innovation_destinations = competitor_innovation_destinations - econparams.n - 1;
            
            firm_state_next(n, :) = (interval_end' == 2).*firm_innovation_destinations + ...
                (interval_end' == 3).*(-competitor_innovation_destinations) + ...
                (interval_end' == 4).*zeros(1, length(firm_strategies)) + ...
                (interval_end' == 5).*firm_state_curr(n, :);
            
        end
        
        firm_state(:,:,i+1) = firm_state_next;
        
    end
        
    out.econparams = econparams;
    out.transition = transition;
    out.sigma_start = sigma_start;
    out.t = t;
    out.firm_profits = squeeze(mean(firm_profits, 1));
            
end

