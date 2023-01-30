%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_transition_stats.m
% Author: Nicholas von Turkovich
% Date: 5/24/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [transition] = compute_transition_stats(econparams, transition, Q0, C0, heta_path, rho_path, final_leg_start, statedep_cutoff)
    
    % Aggregate productivity growth and level
    transition.g = nan(1,size(heta_path,2)+1);
    transition.Q = nan(1,size(heta_path,2)+1);
    transition.Q(1) = Q0;
    
    % Consumption growth and level
    transition.cg = nan(1,size(heta_path,2)+1);
    transition.C = nan(1,size(heta_path,2)+1);
    transition.C(1) = C0;
    
    % Productivity growth
    transition.prodg_outputcomp = nan(1,size(heta_path,2)+1);
    transition.prodg_laborchg = nan(1,size(heta_path,2)+1);
    transition.prodg = nan(1,size(heta_path,2)+1);
    transition.research_productivity = nan(1,size(heta_path,2)+1);
    
    % Gross markup path and sd of markups
    transition.mean_gross_markup = nan(1,size(heta_path,2)+1);
    transition.std_markup = nan(1,size(heta_path,2)+1);
    transition.pct75_gross_markup = nan(1,size(heta_path,2)+1);
    transition.pct90_gross_markup = nan(1,size(heta_path,2)+1);
    
    % Collect labor demand
    transition.L = nan(1,size(heta_path,2)+1);
    transition.L_research = nan(1,size(heta_path,2)+1);
    
    % Welfare vector
    transition.W_path = nan(1,size(heta_path,2)+1);
    transition.muplab = nan(1, size(heta_path,2)+1);
    
    for t = 1:size(heta_path,2)
        
        % Temporary econparams object
        econparams_temp = econparams;

        % Assign object the values for time 't'
        if nargin == 8 && ~isempty(statedep_cutoff)
            econparams_temp = update_heta(econparams_temp, heta_path(:,t), statedep_cutoff);
        else
            econparams_temp = update_heta(econparams_temp, heta_path(t));
        end
        econparams_temp.rho = rho_path(t);
        econparams_temp.value_function = transition.value_function_path(:,t);
        econparams_temp.x_vector = transition.x_vector_path(:,t);
        econparams_temp.x_vector_e = transition.x_vector_e_path(:,t);
        econparams_temp.mu = transition.mu_path(:,t);
        
        % Compute key moments
        econparams_temp = compute_key_moments(econparams_temp, false);
        
        % Productivity growth just depends on current mu
        transition.g(t) = econparams_temp.g;
        transition.Q(t+1) = transition.Q(t)*exp(transition.g(t));
        
        % Gross markup
        transition.mean_gross_markup(t) = econparams_temp.mean_gross_markup;
        transition.pct75_gross_markup(t) = econparams_temp.pct75_gross_markup;
        transition.pct90_gross_markup(t) = econparams_temp.pct90_gross_markup;
        
        % SD of markups
        transition.std_markup(t) = econparams_temp.std_markup;
        
        % Labor demand
        transition.L(t) = econparams_temp.L;
        transition.L_research(t) = econparams_temp.L_research;
        
        % Consumption growth is a function of productivity growth and
        % changes in markups
        transition.cg(t) = transition.g(t) - sum((transition.mu_path(:,t+1) - transition.mu_path(:,t)).*log(econparams.markups((econparams.n + 1):(2*econparams.n + 1))));
        transition.C(t+1) = transition.C(t)*exp(transition.cg(t));
        
        % Compute flow welfare from consumption and labor
        transition.W_path(t) = log(transition.C(t)) - transition.L(t);
        
        % Separate calculation for the distortions associated with markups
        % and labor disutility
        transition.muplab(t) = (-(transition.mu_path(:,t)' * log(econparams.markups(econparams.n+1:2*econparams.n+1))) - transition.L(t));
        
        % Compute productivity growth (i.e., output per unit production
        % labor)
        industry_production_labor = econparams_temp.production_labor_demand((econparams.n + 1):(econparams.n*2 + 1)) + econparams_temp.production_labor_demand((econparams.n + 1):-1:(1));
        
        transition.prodg_outputcomp(t) = sum((transition.mu_path(:,t+1) - transition.mu_path(:,t)).*log(industry_production_labor));
        transition.prodg_laborchg(t) = (sum((transition.mu_path(:,t+1) - transition.mu_path(:,t)).*(industry_production_labor))/...
            sum((transition.mu_path(:,t)).*(industry_production_labor)));
        
        transition.prodg(t) = transition.g(t) + transition.prodg_outputcomp(t) - transition.prodg_laborchg(t);
        
        % Compute aggregate research productivity
        transition.research_productivity(t) = transition.prodg(t)/econparams_temp.L_research;

    end
    
    % Compute welfare including the transition; use BGP formula for last
    % value
    transition.discount = exp(-cumsum([rho_path(1:end)]));
    transition.W = sum(transition.discount(1:length(rho_path)-1).*transition.W_path(1:length(rho_path)-1)) + ...
        transition.discount(length(rho_path))*((log(transition.C(length(rho_path))) - transition.L(length(rho_path)))/rho_path(end) + transition.g(length(rho_path))/(rho_path(end)^2));
    
    % Compute welfare broken down into sub components
    transition.W_prod = sum(transition.discount(1:length(rho_path)-1) .* log(transition.Q(1:length(rho_path)-1))) + transition.discount(length(rho_path))*(log(transition.Q(length(rho_path)))/rho_path(end) + transition.g(length(rho_path))/(rho_path(end)^2));
    transition.W_muplab = sum(transition.discount(1:length(rho_path)-1) .* transition.muplab(1:length(rho_path)-1)) + transition.discount(length(rho_path))*(-(transition.mu_path(:,length(rho_path))' * log(econparams.markups(econparams.n+1:end))) - transition.L(length(rho_path)))/rho_path(end);

    transition.heta_path = heta_path;
    transition.rho_path = rho_path;
    transition.final_leg_start = final_leg_start;

end