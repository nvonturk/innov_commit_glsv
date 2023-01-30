%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_decomp.m
% Author: Nicholas von Turkovich
% Date: 11/19/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FHK] = compute_decomp(econparams_bgp, sim, num_years)

    % Extract simulation variables
    age_firm1 = sim.age_firm1;
    age_firm2 = sim.age_firm2;
    N = sim.N;
    T = sim.T;
    lambda = sim.lambda;
    periods_per_year = sim.periods_per_year;
    innov_firm1 = sim.innov_firm1;
    innov_firm2 = sim.innov_firm2;
    innov_frontier_firm1 = sim.innov_frontier_firm1;
    innov_frontier_firm2 = sim.innov_frontier_firm2;
    state_firm1 = sim.state_firm1;
    state_firm2 = sim.state_firm2;
    
    % Determine starting point for analysis
    T_decomp = T - num_years*periods_per_year + 1;

    % Extract the ages of the firms
    age_firm1_interval = age_firm1(:,(T_decomp:T+1));
    age_firm2_interval = age_firm2(:,(T_decomp:T+1));
    
    % If the firm is of age 0 at any point then entry occurred
    C_firm1 = zeros(size(age_firm1,1),1);
    C_firm2 = C_firm1;
    
    for i = 1:size(C_firm1,1)
        
        C_firm1(i) = ~ismember(0, age_firm1_interval(i,:));
        C_firm2(i) = ~ismember(0, age_firm2_interval(i,:));
        
    end
    
    % Compute the shares of revenue of each firm in the simulation
    if econparams_bgp.kappa == -1
        revenue_share_firm1 = (state_firm1(:,[T_decomp:T+1]) > 0) + ...
            0.5*(state_firm1(:,[T_decomp:T+1]) == 0);
        revenue_share_firm2 = (state_firm2(:,[T_decomp:T+1]) > 0) + ...
            0.5*(state_firm2(:,[T_decomp:T+1]) == 0);
    else
        error("Need to implement for when kappa is finite")
    end
    
    % Compute changes in log productivity for each firm at each point in
    % time
    prod_changes_firm1 = sum(innov_firm1(:,T_decomp:T),2);
    prod_changes_firm2 = sum(innov_firm2(:,T_decomp:T),2);
    
    % Compute sum of changes that pushed the frontier
    prod_changes_frontier_firm1 = sum(innov_frontier_firm1(:,T_decomp:T),2);
    prod_changes_frontier_firm2 = sum(innov_frontier_firm2(:,T_decomp:T),2);
    
    % Compute change in revenue share
    rev_changes_firm1 = revenue_share_firm1(:,end) - revenue_share_firm1(:,1);
    rev_changes_firm2 = revenue_share_firm2(:,end) - revenue_share_firm2(:,1);
    
    WITHIN = zeros(N,1);
    BETWEEN = zeros(N,1);
    CROSS = zeros(N,1);
    ENTRY = zeros(N,1);
    EXIT = zeros(N,1);
    COMP = zeros(N,1);
    TOTAL = zeros(N,1);
    
    for i = 1:N
        
        if C_firm1(i) == 1
            WITHIN(i) = revenue_share_firm1(i,1) * prod_changes_firm1(i);
            BETWEEN(i) = (1 - revenue_share_firm1(i,1)) * (log(lambda) * state_firm1(i,T_decomp)) * (rev_changes_firm1(i));
            CROSS(i) = prod_changes_firm1(i) * rev_changes_firm1(i);
        else
            ENTRY(i) = revenue_share_firm1(i,end)*(prod_changes_firm1(i) + (1 - revenue_share_firm1(i,1))*state_firm1(i,T_decomp)*log(lambda));
        end
        
        if C_firm2(i) == 1
            WITHIN(i) = WITHIN(i) + revenue_share_firm2(i,1) * prod_changes_firm2(i);
            BETWEEN(i) = BETWEEN(i) + (1 - revenue_share_firm2(i,1)) * (log(lambda) * state_firm2(i,T_decomp)) * (rev_changes_firm2(i));
            CROSS(i) = CROSS(i) + prod_changes_firm2(i) * rev_changes_firm2(i);
        else
            ENTRY(i) = ENTRY(i) + revenue_share_firm2(i,end)*(prod_changes_firm2(i) + (1 - revenue_share_firm2(i,1))*state_firm2(i,T_decomp)*log(lambda));
        end
        
        COMP(i) = prod_changes_frontier_firm1(i) + prod_changes_frontier_firm2(i);
        TOTAL(i) = WITHIN(i) + BETWEEN(i) + CROSS(i) + ENTRY(i);
        
        if abs(TOTAL(i) - COMP(i)) > 1e-10
            error("FHK broken")
        end
        
    end
    
    FHK.WITHIN = WITHIN;
    FHK.BETWEEN = BETWEEN;
    FHK.CROSS = CROSS;
    FHK.ENTRY = ENTRY;
    FHK.EXIT = EXIT;  
    
    % See equation 2
    
%     % Growth in unweighted average productivity level
%     growth_unweighted_average = (C_firm1 .* prod_changes_firm1 + ...
%         C_firm2 .* prod_changes_firm2) ./ (C_firm1 + C_firm2);
%     
%     % If NaN then divided by 0 -> no continuing firms
%     growth_unweighted_average(isnan(growth_unweighted_average)) = 0;
%     
%     % Covariance change
%     
%     % Compute covariance terms for each firm and time interval
%     cov_before_firm1 = (revenue_share_firm1(:,1) - 0.5) .* (0.5*log(lambda)*state_firm1(:,T_decomp));
%     cov_before_firm2 = (revenue_share_firm2(:,1) - 0.5) .* (0.5*log(lambda)*state_firm2(:,T_decomp));
%     cov_after_firm1 = (revenue_share_firm1(:,end) - 0.5) .* (0.5*log(lambda)*state_firm1(:,T+1));
%     cov_after_firm2 = (revenue_share_firm2(:,end) - 0.5) .* (0.5*log(lambda)*state_firm2(:,T+1));
%     
%     cov_change = cov_after_firm1 + cov_after_firm2 - ...
%         (cov_before_firm1 + cov_before_firm2);
%     
%     % Since if there is an entry, the "continuing set" is just a single
%     % firm so their log productivity level will be the mean for the
%     % continuing set; all covariance terms are 0
%     cov_change = (C_firm1 & C_firm2) .* cov_change;
%     
%     % Entry term
%     dopd_entry = (~C_firm1 & C_firm2).*(revenue_share_firm1(:,end)).*(log(lambda) .* state_firm1(:,T+1)) + ...
%         (~C_firm2 & C_firm1).*(revenue_share_firm2(:,end)).*(log(lambda) .* state_firm2(:,T+1));
%     
%     % Exit term
%     dopd_exit = (~C_firm1 & C_firm2).*(revenue_share_firm1(:,1)).*(log(lambda) .* state_firm2(:,T_decomp)) + ...
%         (~C_firm2 & C_firm1).*(revenue_share_firm2(:,1)).*(log(lambda) .* state_firm1(:,T_decomp));
%     
%     % Entry add on when both firms are entrants
%     dopd_entry = dopd_entry + (~C_firm1 & ~C_firm2).*(prod_changes_frontier_firm1 + prod_changes_frontier_firm2);
%     
%     if sum(abs(growth_unweighted_average + cov_change + dopd_entry + dopd_exit - COMP), 'all') > 1e-10
%         error("OP decomposition broken")
%     end
%     
%     OP.WITHIN = growth_unweighted_average;
%     OP.COV = cov_change;
%     OP.NET_ENTRY = dopd_entry + dopd_exit;
%     
end