%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: simulate_economy.m
% Author: Nicholas von Turkovich
% Date: 11/8/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams_bgp] = simulate_economy(econparams_bgp, flag_check)

    %% Extract Variables for Loops to Cut Down on Runtime
    
    % Model variables
    x_vector = econparams_bgp.x_vector;
    x_vector_e = econparams_bgp.x_vector_e;
    heta_vector = econparams_bgp.heta_vector;
    mu = econparams_bgp.mu;
    innovation_transitions = econparams_bgp.innovation_transitions;
    innovation_transitions_e = econparams_bgp.innovation_transitions_e;
    patent_transitions = econparams_bgp.patent_transitions;
    lambda = econparams_bgp.lambda;
    entry_flag = econparams_bgp.entry_flag;
    total_labor_demand = econparams_bgp.total_labor_demand;
    value_function = econparams_bgp.value_function;
    
    %% Setup
        
    % Number of firms
    N = 5*10^4;
    
    % Determine if the time scale needs to be adjusted
    if (2*max(x_vector) + max(x_vector_e) + max(heta_vector)) > 1
        scale_time = 0.99/(2*max(x_vector) + max(x_vector_e) + max(heta_vector));
    else
        scale_time = 1;
    end
    
    % Periods per year (larger than 1/frequency if parameterization leads
    % to more than 1 event in expectation each time period)
    periods_per_year = ceil(1/econparams_bgp.frequency/scale_time);
    
    % Adjust rates for scale time
    x_vector = x_vector*scale_time;
    x_vector_e = x_vector_e*scale_time;
    heta_vector = heta_vector*scale_time;
    
    % Xi parameter (this is the gross rate of technology atrophy each
    % period)
    xi = exp(log(1)/periods_per_year);
    
    % Number of discrete intervals
    T = 12*periods_per_year;
    
    % Generate draw of events
    seed = 1000;
    rng(seed, 'Threefry')
    draws = rand(N, T);
    
    % Data structures for tracking state of industries and firms
    
    % Tracking state of industry and firms (T+1 since including starting state)
    state_industry = nan(N, T+1);
    state_firm1 = nan(N, T+1);
    state_firm2 = nan(N, T+1);
    
    % Tracking innovation changes at each time step of either firm (also
    % includes entrants that take over a half line); frontier innovation
    % only if it pushes the productivity level of the industry
    innov_firm1 = ones(N, T)*log(xi);
    innov_firm2 = ones(N, T)*log(xi);
    innov_frontier_firm1 = zeros(N, T);
    innov_frontier_firm2 = zeros(N, T);
    
    % Tracks changes in value from innovations for matching Kogan et al.
    % moments
    val_changes_firm1 = zeros(N, T);
    val_changes_firm2 = zeros(N, T);
    
    % Track which firm is the leader (T+1 since including starting state)
    leader_tracker = [ones(N,1), nan(N, T)];
    
    % Track what type of event occurred (none, leader innovation, 
    % laggard innovation, entry, patent expiry)
    event_tracker = zeros(N,T);

    % Tracking age of firms (start at Inf and drop to 0 when there is
    % entry; starting at T+1 since including starting state)
    age_firm1 = Inf(N,T+1);
    age_firm2 = Inf(N,T+1);
    
    % Instantiate set of firms and industries
    
    % Starting state of economy (take distribution and compute number of
    % industries per gap
    init_industry_gaps = nan(N,1);
    num_industries = round(abs(mu*N));

    % Loop through and give the industry gap for each block of N firms
    iter = 1;
    for s = 1:length(mu)
       
        init_industry_gaps(iter:(num_industries(s) + iter - 1)) = s-1;
        iter = iter + num_industries(s);
        
    end
    
    % Check to make sure that rounded amount is closest to original number
    % of firms
    if abs(sum(num_industries) - N)/N > 0.001
        error("Rounding introduced noticeable differences")
    end
    
    % If rounding introduces difference between N and number of firms
    % calculated, fill in / remove remaining industries with gap associated with
    % mode of distribution
    if sum(num_industries) < N
                
        s_to_add_mass = find(mu == max(mu)) - 1;
        init_industry_gaps(isnan(init_industry_gaps)) = s_to_add_mass;
        
    elseif sum(num_industries) > N
                
        s_to_remove_mass = find(mu == max(mu)) - 1;
        indices_to_remove = find(init_industry_gaps == s_to_remove_mass, sum(num_industries) - N);
        init_industry_gaps(indices_to_remove) = [];
        
    end
    
    % Initialize states
    state_industry(:,1) = init_industry_gaps;
    state_firm1(:,1) = init_industry_gaps;
    state_firm2(:,1) = -init_industry_gaps;
    
    %% Simulation of Industry Events
    
    % For each time period and for each industry, 
    for t = 1:T
               
        for n = 1:N
            
            % Determine state of industry n at time t
            state = state_industry(n, t);
            
            % Extract innovation and patent expiry rates
            x_leader = x_vector(econparams_bgp.n + 1 + state);
            x_follower = x_vector(econparams_bgp.n + 1 - state);
            heta = heta_vector(econparams_bgp.n + 1 + state);
            
            % If directed entry, extract the appropriate entry innovation
            % rate
            if entry_flag == -1
                x_entrant = x_vector_e;
            else
                x_entrant = x_vector_e(econparams_bgp.n + 1 - state);
            end
            
            % Compute probability of no event
            no_event = 1 - x_leader - x_follower - heta - x_entrant;
            
            if no_event < 0 || no_event > 1
               error("Negative event or increase frequency"); 
            end
            
            % Determine what happens given uniform prob. draw
            draw = draws(n, t);
            
            if draw <= no_event
                
                % Industry stays where it is
                
                state_industry(n, t+1) = state;
                state_firm1(n, t+1) = state_firm1(n, t);
                state_firm2(n, t+1) = state_firm2(n, t);
                leader_tracker(n, t+1) = leader_tracker(n, t);
                age_firm1(n, t+1) = age_firm1(n, t)+1;
                age_firm2(n,t+1) = age_firm2(n,t)+1;
                
                if leader_tracker(n,t) == 1
                    
                    innov_frontier_firm1(n,t) = log(xi);
                    
                else
                    
                    innov_frontier_firm2(n,t) = log(xi);
                    
                end
                
            elseif draw > no_event && draw <= (no_event + x_leader)
                
                % Leader has innovated
                event_tracker(n,t) = 1;
                
                % Determine probability distribution over new states
                pdist = innovation_transitions(:, econparams_bgp.n + 1 + state);
                
                % Set destination state
                new_state = NaN;
                
                % Loop through states to find out how much the leader
                % innovated
                interval_start = no_event;
                
                for i = 1:length(pdist)
                    
                    interval_end = interval_start + pdist(i)*x_leader;
                    
                    if draw > interval_start && draw <= interval_end
                        
                        if ~isnan(new_state)
                           error("You already set the state"); 
                        end
                        
                        % This is where the leader is headed
                        new_state = abs(i - econparams_bgp.n - 1);
                        break
                        
                    end
                    
                    interval_start = interval_end;
                    
                end
                
                if isnan(new_state)
                    save("../../debug.mat")
                    error("Need to set the state")
                end
                
                state_industry(n,t+1) = new_state;
                
                % need to determine which of the two firms to advance
                if (leader_tracker(n, t) == 1)
                                        
                    state_firm1(n, t+1) = i - econparams_bgp.n - 1;
                    state_firm2(n, t+1) = -(i - econparams_bgp.n - 1);
                    
                    innov_firm1(n,t) = innov_firm1(n,t) + (state_firm1(n, t+1) - state_firm1(n, t))*log(lambda);
                    innov_frontier_firm1(n,t) = innov_firm1(n,t);
                    
                    val_changes_firm1(n, t) = value_function(econparams_bgp.n + 1 + state_firm1(n, t+1)) - value_function(econparams_bgp.n + 1 + state_firm1(n, t));
                    
                    leader_tracker(n, t+1) = 1;
                    
                else
                    
                    state_firm2(n,t+1) = i - econparams_bgp.n - 1;
                    state_firm1(n,t+1) = -(i - econparams_bgp.n - 1);
                    
                    innov_firm2(n,t) = innov_firm2(n,t) + (state_firm2(n, t+1) - state_firm2(n, t))*log(lambda);
                    innov_frontier_firm2(n,t) = innov_firm2(n,t);
                    
                    val_changes_firm2(n, t) = value_function(econparams_bgp.n + 1 + state_firm2(n, t+1)) - value_function(econparams_bgp.n + 1 + state_firm2(n, t));
                    
                    leader_tracker(n, t+1) = 2;
                    
                end
                
                age_firm1(n,t+1) = age_firm1(n,t)+1;
                age_firm2(n,t+1) = age_firm2(n,t)+1;
                
            elseif draw > (no_event + x_leader) && draw <= (no_event + x_leader + x_follower)
                
                % follower has innovated
                event_tracker(n,t) = 2;
                
                % determine probability distribution over new states
                pdist = innovation_transitions(:, econparams_bgp.n + 1 - state);
                
                % set destination state
                new_state = NaN;
                
                % loop through states to find out how much the leader
                % innovated
                interval_start = no_event + x_leader;
                
                for i = 1:length(pdist)
                    
                    interval_end = interval_start + pdist(i)*x_follower;
                    
                    if draw > interval_start && draw <= interval_end
                        
                        if ~isnan(new_state)
                           error("You already set the state"); 
                        end
                        
                        % this is where the follower is headed
                        new_state = abs(i - econparams_bgp.n - 1);
                        break
                        
                    end
                    
                    interval_start = interval_end;
                    
                end
                
                if isnan(new_state)
                    save("../../debug.mat")
                    error("Need to set the state")
                end
                
                state_industry(n,t+1) = new_state;
                
                if (leader_tracker(n, t) == 1)
                                        
                    state_firm2(n, t+1) = i - econparams_bgp.n - 1;
                    state_firm1(n, t+1) = -(i - econparams_bgp.n - 1);
                    
                    innov_firm2(n, t) = innov_firm2(n,t) + (state_firm2(n, t+1) - state_firm2(n,t))*log(lambda); 
                    innov_frontier_firm2(n, t) = max(state_firm2(n, t+1)*log(lambda),0); 
                    
                    val_changes_firm2(n, t) = value_function(econparams_bgp.n + 1 + state_firm2(n, t+1)) - value_function(econparams_bgp.n + 1 + state_firm2(n, t));

                    if i > econparams_bgp.n + 1
                        leader_tracker(n, t+1) = 2;
                        innov_frontier_firm2(n,t) = innov_frontier_firm2(n,t) + log(xi);
                    else
                        leader_tracker(n, t+1) = 1;
                        innov_frontier_firm1(n,t) = innov_frontier_firm1(n,t) + log(xi);
                    end
                                        
                else
                                        
                    state_firm1(n,t+1) = i - econparams_bgp.n - 1;
                    state_firm2(n,t+1) = -(i - econparams_bgp.n - 1);
                    
                    innov_firm1(n, t) = innov_firm1(n,t) + (state_firm1(n, t+1) - state_firm1(n,t))*log(lambda); 
                    innov_frontier_firm1(n, t) = max(state_firm1(n, t+1)*log(lambda), 0); 
                    
                    val_changes_firm1(n, t) = value_function(econparams_bgp.n + 1 + state_firm1(n, t+1)) - value_function(econparams_bgp.n + 1 + state_firm1(n, t));

                    if i > econparams_bgp.n + 1
                        leader_tracker(n,t+1) = 1;
                        innov_frontier_firm1(n,t) = innov_frontier_firm1(n,t) + log(xi);
                    else
                        leader_tracker(n,t+1) = 2;
                        innov_frontier_firm2(n,t) = innov_frontier_firm2(n,t) + log(xi);
                    end
                                        
                end
                
                age_firm1(n,t+1) = age_firm1(n,t)+1;
                age_firm2(n,t+1) = age_firm2(n,t)+1;
                
                
            elseif draw > (no_event + x_leader + x_follower) && draw <= (no_event + x_leader + x_follower + x_entrant)
                
                % entrant has innovated
                event_tracker(n,t) = 3;
                
                % determine probability distribution over new states
                pdist = innovation_transitions_e(:, econparams_bgp.n + 1 - state);  
                
                % set destination state
                new_state = NaN;
                
                % loop through states to find out how much the leader
                % innovated
                interval_start = no_event + x_leader + x_follower;
                
                for i = 1:length(pdist)
                    
                    interval_end = interval_start + pdist(i)*x_entrant;
                    
                    if draw > interval_start && draw <= interval_end
                        
                        if ~isnan(new_state)
                           error("You already set the state"); 
                        end
                        
                        % this is where the follower is headed
                        new_state = abs(i - econparams_bgp.n - 1);
                        break
                        
                    end
                    
                    interval_start = interval_end;
                    
                end
                
                if isnan(new_state)
                    save("../../debug.mat")
                    error("Need to set the state")
                end
                
                state_industry(n,t+1) = new_state;
                
                if state == 0
                    
                    if draw > (no_event + x_leader + x_follower + x_entrant/2)
                       
                        % knock out firm 2
                        state_firm2(n,t+1) = i - econparams_bgp.n - 1;
                        state_firm1(n,t+1) = -(i - econparams_bgp.n - 1);
                        innov_firm2(n,t) = innov_firm2(n,t) + state_firm2(n,t+1)*log(lambda);
                        innov_frontier_firm2(n,t) = innov_firm2(n,t);
                        
                        val_changes_firm2(n, t) = value_function(econparams_bgp.n + 1 + state_firm2(n, t+1)) - value_function(econparams_bgp.n + 1 + state_firm2(n, t));
                        val_changes_firm2(n,1:(t-1)) = 0;
                        
                        leader_tracker(n,t+1) = 2;
                        age_firm1(n,t+1) = age_firm1(n,t)+1;
                        age_firm2(n,t+1) = 0;
                        
                        
                    else
                        
                        % knock out firm 1
                        state_firm1(n,t+1) = i - econparams_bgp.n - 1;
                        state_firm2(n,t+1) = -(i - econparams_bgp.n - 1);
                        innov_firm1(n,t) = innov_firm1(n,t) + state_firm1(n,t+1)*log(lambda);
                        innov_frontier_firm1(n,t) = innov_firm1(n,t);
                        
                        val_changes_firm1(n, t) = value_function(econparams_bgp.n + 1 + state_firm1(n, t+1)) - value_function(econparams_bgp.n + 1 + state_firm1(n, t));
                        val_changes_firm1(n,1:(t-1)) = 0;
                        
                        leader_tracker(n,t+1) = 1;
                        age_firm1(n,t+1) = 0;
                        age_firm2(n,t+1) = age_firm2(n,t)+1;
                        
                    end
                    
                    
                    
                elseif state ~= 0 && (leader_tracker(n, t) == 1)
                    
                    state_firm2(n, t+1) = i - econparams_bgp.n - 1;
                    state_firm1(n, t+1) = -(i - econparams_bgp.n - 1);
                    innov_firm2(n,t) = innov_firm2(n,t) + (state_firm2(n,t+1) - state_firm2(n,t))*log(lambda);
                    innov_frontier_firm2(n,t) = max(state_firm2(n,t+1)*log(lambda),0);
                    
                    val_changes_firm2(n, t) = value_function(econparams_bgp.n + 1 + state_firm2(n, t+1)) - value_function(econparams_bgp.n + 1 + state_firm2(n, t));
                    val_changes_firm2(n,1:(t-1)) = 0;
                    
                    if i > econparams_bgp.n + 1
                        leader_tracker(n, t+1) = 2;
                        innov_frontier_firm2(n,t) = innov_frontier_firm2(n,t) + log(xi);
                    else
                        leader_tracker(n, t+1) = 1;
                        innov_frontier_firm1(n,t) = innov_frontier_firm1(n,t) + log(xi);
                    end
                    
                    age_firm1(n,t+1) = age_firm1(n,t)+1;
                    age_firm2(n,t+1) = 0;
                                        
                elseif state ~= 0 && (leader_tracker(n, t) == 2)
                    
                    state_firm1(n,t+1) = i - econparams_bgp.n - 1;
                    state_firm2(n,t+1) = -(i - econparams_bgp.n - 1);
                    innov_firm1(n,t) = innov_firm1(n,t) + (state_firm1(n,t+1) - state_firm1(n,t))*log(lambda);
                    innov_frontier_firm1(n,t) = max(state_firm1(n,t+1)*log(lambda),0);
                    
                    val_changes_firm1(n, t) = value_function(econparams_bgp.n + 1 + state_firm1(n, t+1)) - value_function(econparams_bgp.n + 1 + state_firm1(n, t));
                    val_changes_firm1(n,1:(t-1)) = 0;

                    if i > econparams_bgp.n + 1
                        leader_tracker(n,t+1) = 1;
                        innov_frontier_firm1(n,t) = innov_frontier_firm1(n,t) + log(xi);
                    else
                        leader_tracker(n,t+1) = 2;
                        innov_frontier_firm2(n,t) = innov_frontier_firm2(n,t) + log(xi);
                    end
                    
                    age_firm1(n,t+1) = 0;
                    age_firm2(n,t+1) = age_firm2(n,t)+1;
                                        
                end
                
                                
            elseif draw > (no_event + x_leader + x_follower + x_entrant)
            
                % patent expired
                event_tracker(n,t) = 4;
                
                if state == 0
                    error("Shouldn't have PE in neck and neck")
                end
                
                % determine probability distribution over new states for
                % leader
                pdist = patent_transitions(:, econparams_bgp.n + 1 + state);
                
                % set destination state
                new_state = NaN;
                
                interval_start = no_event + x_leader + x_follower + x_entrant;
                
                for i = 1:length(pdist)
                    
                    interval_end = interval_start + pdist(i)*heta;
                    
                    if draw > interval_start && draw <= interval_end
                        
                        if ~isnan(new_state)
                           error("You already set the state"); 
                        end
                        
                        % this is where the follower is headed
                        new_state = abs(i - econparams_bgp.n - 1);
                        break
                        
                    end
                    
                    interval_start = interval_end;
                    
                end
                
                if isnan(new_state)
                    save("../../debug.mat")
                    error("Need to set the state")
                end
                
                state_industry(n,t+1) = new_state;
                
                % since from perspective of leader, leader drops to i and
                % the follower drops to complement of i
                
                if (leader_tracker(n, t) == 1)
                    
                    state_firm1(n, t+1) = i - econparams_bgp.n - 1;
                    state_firm2(n, t+1) = -(i - econparams_bgp.n - 1);
                    innov_firm2(n, t) = innov_firm2(n,t) + (state_firm2(n,t+1) - state_firm2(n,t))*log(lambda);
                    
                    if i < econparams_bgp.n + 1
                        leader_tracker(n,t+1) = 2;
                        innov_frontier_firm2(n, t) = log(xi);
                    else
                        leader_tracker(n,t+1) = 1;
                        innov_frontier_firm1(n, t) = log(xi);
                    end
                    
                else
                    
                    state_firm2(n, t+1) = i - econparams_bgp.n - 1;
                    state_firm1(n, t+1) = -(i - econparams_bgp.n - 1);
                    innov_firm1(n, t) = innov_firm1(n,t) + (state_firm1(n,t+1) - state_firm1(n,t))*log(lambda);
                    
                    if i < econparams_bgp.n + 1
                        leader_tracker(n,t+1) = 1;
                        innov_frontier_firm1(n, t) = log(xi);
                    else
                        leader_tracker(n,t+1) = 2;
                        innov_frontier_firm2(n, t) = log(xi);
                    end
                    
                end
                
                age_firm1(n,t+1) = age_firm1(n,t)+1;
                age_firm2(n,t+1) = age_firm2(n,t)+1;
                
                
            end
        
        end

    end
    
    %% Checks on Simulation
    
    if flag_check
    
        % leader and follower position should sum to 2*s_max + 2
        check1 = state_firm1 + state_firm2;
        if ~isempty(check1(check1 ~= 0))
           error("Wrong states") 
        end

        % check that state_firm1 and state_firm2 distance to state is the same
        if sum(abs(state_firm1) - state_industry, 'all') ~= 0
           error("Wrong states 2") 
        end

        % leader is always at least in position 51
        check3 = (leader_tracker == 1).*state_firm1 + (leader_tracker == 2).*state_firm2;

        if min(check3(:)) < 0 || max(check3(:)) > econparams_bgp.n
            error("Wrong states 3")
        end

        % check to see if check3 is always equal to state_industry
        if ~isempty(state_industry(state_industry ~= check3))
            error("Wrong states 4")
        end
        
        % Check that innovation changes match state changes
        if sum(abs(innov_firm1 - innov_firm2) ~= abs(state_firm1(:,2:end) - state_firm1(:,1:(end-1)))*log(lambda), 'all')
            error("Wrong innovation counting 1")
        end

        state_space_start_finish = [state_industry(:,1), state_industry(:,round(T/2)), state_industry(:,end)];
        age_start_finish = [[age_firm1(:,1), age_firm1(:,round(T/2)), age_firm1(:,round(3*T/4)), age_firm1(:,end)];...
            [age_firm2(:,1), age_firm2(:,round(T/2)), age_firm2(:,round(3*T/4)), age_firm2(:,end)]];

        plot_global_settings();

        f = figure('name', 'Simulation Check 1');
        hist(state_space_start_finish, max(state_space_start_finish(:)) + 1)
        xlabel("$\sigma$")
        ylabel("Count")
        legend("Beginning of Simulation", "Middle of Simulation", "End of Simulation")
        saveas(f,"../Output/Figures/Robustness/" + f.Name, 'eps');

        f = figure('name', 'Simulation Check 2')
        histogram(state_space_start_finish(:,3), 'BinMethod', 'integers', 'Normalization', 'probability')
        ylim([0, 0.5])
        hold on
        yyaxis right
        ylim([0, 0.5])
        plot(0:econparams_bgp.n, mu, 'k-')
        hold off
        yyaxis left
        xlabel("$\sigma$")
        ylabel("Density")
        legend("Simulation", "BGP Density");
        saveas(f,"../Output/Figures/Robustness/" + f.Name, 'eps');
        
        f = figure('name', 'Simulation Check 3 - Entry Age for Firms 1,2')
        hist([age_firm1(:,T), age_firm2(:,T)]*econparams_bgp.frequency)
        xlabel("Age (Years)")
        ylabel("Density")
        legend("Firm 1", "Firm 2")
        saveas(f,"../Output/Figures/Robustness/" + f.Name, 'eps');
        
        f = figure('name', 'Simulation Check 4 - Distribution of Ages')
        hist(age_start_finish*econparams_bgp.frequency, 20)
        xlabel("Age (Years)")
        ylabel("Count")
        legend("Start of Simulation", "T/2", "3T/4", "T")
        saveas(f,"../Output/Figures/Robustness/" + f.Name, 'eps');
        
    end
    
    %% Gather Simulation Objects
    sim.age_firm1 = age_firm1;
    sim.age_firm2 = age_firm2;
    sim.N = N;
    sim.T = T;
    sim.lambda = lambda;
    sim.periods_per_year = periods_per_year;
    sim.innov_firm1 = innov_firm1;
    sim.innov_firm2 = innov_firm2;
    sim.innov_frontier_firm1 = innov_frontier_firm1;
    sim.innov_frontier_firm2 = innov_frontier_firm2;
    sim.state_firm1 = state_firm1;
    sim.state_firm2 = state_firm2;
    sim.state_industry = state_industry;
    sim.leader_tracker = leader_tracker;
    sim.event_tracker = event_tracker;
    
    %% FHK 
    
    [FHK] = compute_decomp(econparams_bgp, sim, 5);
    
    econparams_bgp.FHK_WITHIN_adj_sh_5y = sum(FHK.WITHIN)/(sum(FHK.WITHIN) + sum(FHK.BETWEEN) + sum(FHK.CROSS));
    econparams_bgp.FHK_WITHIN_sh_5y = sum(FHK.WITHIN)/(sum(FHK.WITHIN) + sum(FHK.BETWEEN) + sum(FHK.CROSS) + sum(FHK.ENTRY) + sum(FHK.EXIT));
    econparams_bgp.FHK_BETWEEN_sh_5y = sum(FHK.BETWEEN)/(sum(FHK.WITHIN) + sum(FHK.BETWEEN) + sum(FHK.CROSS) + sum(FHK.ENTRY) + sum(FHK.EXIT));
    econparams_bgp.FHK_CROSS_sh_5y = sum(FHK.CROSS)/(sum(FHK.WITHIN) + sum(FHK.BETWEEN) + sum(FHK.CROSS) + sum(FHK.ENTRY) + sum(FHK.EXIT));
    econparams_bgp.FHK_ENTRY_sh_5y = sum(FHK.ENTRY)/(sum(FHK.WITHIN) + sum(FHK.BETWEEN) + sum(FHK.CROSS) + sum(FHK.ENTRY) + sum(FHK.EXIT));
    econparams_bgp.FHK_EXIT_sh_5y = sum(FHK.EXIT)/(sum(FHK.WITHIN) + sum(FHK.BETWEEN) + sum(FHK.CROSS) + sum(FHK.ENTRY) + sum(FHK.EXIT));
    
    econparams_bgp.FHK_WITHIN_pp_5y = mean(FHK.WITHIN);
    econparams_bgp.FHK_BETWEEN_pp_5y = mean(FHK.BETWEEN);
    econparams_bgp.FHK_CROSS_pp_5y = mean(FHK.CROSS);
    econparams_bgp.FHK_ENTRY_pp_5y = mean(FHK.ENTRY);    
    econparams_bgp.FHK_EXIT_pp_5y = mean(FHK.EXIT);

    %% Employment Share

    t_5 = round(6/econparams_bgp.frequency/scale_time + 1);
    
    young5_firm1 = age_firm1(:,t_5) < Inf;
    young5_firm2 = age_firm2(:,t_5) < Inf;
    
    total5_employment_firm1 = total_labor_demand(state_firm1(:,t_5) + econparams_bgp.n + 1);
    total5_employment_firm2 = total_labor_demand(state_firm2(:,t_5) + econparams_bgp.n + 1);
    
    econparams_bgp.emp5 = (sum(total5_employment_firm1(young5_firm1)) + sum(total5_employment_firm2(young5_firm2)))/(sum(total5_employment_firm1) + sum(total5_employment_firm2));
    
    t_10 = round(11/econparams_bgp.frequency/scale_time + 1);
    
    young10_firm1 = age_firm1(:,t_10) < Inf;
    young10_firm2 = age_firm2(:,t_10) < Inf;
    
    total10_employment_firm1 = total_labor_demand(state_firm1(:,t_10) + econparams_bgp.n + 1);
    total10_employment_firm2 = total_labor_demand(state_firm2(:,t_10) + econparams_bgp.n + 1);
    
    econparams_bgp.emp10 = (sum(total10_employment_firm1(young10_firm1)) + sum(total10_employment_firm2(young10_firm2)))/(sum(total10_employment_firm1) + sum(total10_employment_firm2));

    %% Innovation Output
    
    % Cumulate value changes over designated time span
    pct_value_change_firm1 = zeros(N,T/periods_per_year);
    pct_value_change_firm2 = zeros(N,T/periods_per_year);
    
    for i = 1:size(pct_value_change_firm1,2)
        
        eoy = i*periods_per_year + 1;
        pct_value_change_firm1(:,i) = sum(val_changes_firm1(:,(eoy - periods_per_year):(eoy - 1)),2) ./ value_function(econparams_bgp.n + 1 + state_firm1(:,eoy));
        pct_value_change_firm2(:,i) = sum(val_changes_firm2(:,(eoy - periods_per_year):(eoy - 1)),2) ./ value_function(econparams_bgp.n + 1 + state_firm2(:,eoy));
        
    end
    
    pct_value_change_total = [pct_value_change_firm1; pct_value_change_firm2];
    
    econparams_bgp.innovation_output_mean = mean(pct_value_change_total(:));
    econparams_bgp.innovation_output_median = prctile(pct_value_change_total(:), 50);
    econparams_bgp.innovation_output_pct90 = prctile(pct_value_change_total(:), 90);
    
    %% R&D to Sales
    
    % Calculate sales for each firm at each time period
    if ~econparams_bgp.entry_flag
    
        if econparams_bgp.kappa == -1
            sales_vec = [repmat(0,1,econparams_bgp.n), 0.5, repmat(1,1,econparams_bgp.n)];
        else
            error("Not set up for kappa ~= -1")
        end
        
        % Calculate the profit of each firm at each time period
        profit_firm1 = econparams_bgp.profits(state_firm1 + econparams_bgp.n + 1);
        profit_firm2 = econparams_bgp.profits(state_firm2 + econparams_bgp.n + 1);
        
        % Calculate the sales of each firm at each time period
        sales_firm1 = sales_vec(state_firm1 + econparams_bgp.n + 1);
        sales_firm2 = sales_vec(state_firm2 + econparams_bgp.n + 1);

        % Calculate the R&D of each firm at each time period
        rd_firm1 = econparams_bgp.research_labor_demand(state_firm1 + econparams_bgp.n + 1)*econparams_bgp.omega;
        rd_firm2 = econparams_bgp.research_labor_demand(state_firm2 + econparams_bgp.n + 1)*econparams_bgp.omega;
        
        % Reshape to annualize the data
        profit_firm1_ann = reshape(profit_firm1(:,1:T), size(profit_firm1,1), periods_per_year, T/periods_per_year);
        profit_firm2_ann = reshape(profit_firm2(:,1:T), size(profit_firm2,1), periods_per_year, T/periods_per_year);
        sales_firm1_ann = reshape(sales_firm1(:,1:T), size(sales_firm1,1), periods_per_year, T/periods_per_year);
        sales_firm2_ann = reshape(sales_firm2(:,1:T), size(sales_firm2,1), periods_per_year, T/periods_per_year);
        rd_firm1_ann = reshape(rd_firm1(:,1:T), size(rd_firm1,1), periods_per_year, T/periods_per_year);
        rd_firm2_ann = reshape(rd_firm2(:,1:T), size(rd_firm2,1), periods_per_year, T/periods_per_year);
        
        % Sum over months in a year and then reshape again
        profit_firm1_ann = reshape(nansum(profit_firm1_ann,2),size(profit_firm1_ann,1),size(profit_firm1_ann,3));
        profit_firm2_ann = reshape(nansum(profit_firm2_ann,2),size(profit_firm2_ann,1),size(profit_firm2_ann,3));
        sales_firm1_ann = reshape(nansum(sales_firm1_ann,2),size(sales_firm1_ann,1),size(sales_firm1_ann,3));
        sales_firm2_ann = reshape(nansum(sales_firm2_ann,2),size(sales_firm2_ann,1),size(sales_firm2_ann,3));
        rd_firm1_ann = reshape(nansum(rd_firm1_ann,2),size(rd_firm1_ann,1),size(rd_firm1_ann,3));
        rd_firm2_ann = reshape(nansum(rd_firm2_ann,2),size(rd_firm2_ann,1),size(rd_firm2_ann,3));
        
        % Profits for all firms
        profits_ann = [profit_firm1_ann;...
            profit_firm2_ann];
        
        % Compute R&D to Sales
        rd_sales = [rd_firm1_ann./sales_firm1_ann;...
            rd_firm2_ann./sales_firm2_ann];
        
        % If entry is Inf, then place a nan
        rd_sales(isinf(rd_sales)) = nan;
        
        % Unconditional mean and median
        econparams_bgp.rd_sales_mean = nanmean(nanmean(rd_sales,1));
        econparams_bgp.rd_sales_median = nanmean(nanmedian(rd_sales,1));
        
        % Quintile median
        rd_sales_quintile_median = nan(5,size(profits_ann,2));
        
        % Break down by quintiles for annual profit
        for (i = 1:size(profits_ann,2))
            
            % Identify firms with valid RD_sales/profits
            idx = find(profits_ann(:,i) >= 0 & ...
                ~isnan(rd_sales(:,i)));
            
            % Extract values for those firms
            rd_sales_iter = rd_sales(idx,i);
            profits_ann_iter = profits_ann(idx,i);
            
            % Sort the firms in ascending order by profit
            [~, sort_idx] = sort(profits_ann_iter, 'ascend');
            
            % Calculate break point indices for quintiles
            segment_length = ceil(length(sort_idx)/5);
            idx_iter = 1;
            for j = 1:5
                
                idx_iter_end = min(idx_iter + segment_length,length(sort_idx));
                
                rd_sales_quintile_median(j,i) = median(rd_sales_iter(sort_idx(idx_iter:idx_iter_end)));
                
                idx_iter = idx_iter_end + 1;
                
            end
            
        end
        
        % Extract R&D to sales for only those in the top profit bracket
        econparams_bgp.rd_sales_topprof_median = mean(rd_sales_quintile_median(end,:));
        

        
    end
    
    
            
end

