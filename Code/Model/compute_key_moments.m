%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_key_moments.m
% Author: Nicholas von Turkovich
% Date: 10/25/2021
% Note(s): Produces moments of interest given a solved model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = compute_key_moments(econparams, fhk_flag)
        
    % Compute growth rate
    econparams = compute_growth_rate(econparams);
    
    % Compute mean gross markup
    econparams.mean_gross_markup = sum(econparams.mu .* econparams.markups((econparams.n + 1):(2*econparams.n + 1)));
    
    % Compute standard deviation of markups
    econparams.std_markup = std(econparams.markups((econparams.n + 1):(2*econparams.n + 1)), abs(econparams.mu));
    
    % Gross percentiles for markups
    
    % Interplation switch
    interp_switch = 1;
    
    if interp_switch
    
        % Find the industry gap where the density reaches 100 (hence adjustment
        % by -1)
        reaches_100 = find(cumsum(econparams.mu) >= 1 - 1e-10, 1) - 1;
        if isempty(reaches_100)
            error("Why does the sum of mu not equal 1?")
        end
        
        % Find the first industry gap where there is positive density
        above_0 = find(cumsum(econparams.mu) >= 1e-10, 1) - 1;
        if isempty(above_0)
            error("No positive entry in mu")
        end
        
        if above_0 == reaches_100
            % All of mass is concentrated in a single industry gap
            gap_50 = above_0;
            gap_75 = above_0;
            gap_90 = above_0;
            gap_95 = above_0;
            gap_99 = above_0;
        elseif reaches_100 > above_0
            
            % There is at least two gaps with density however if the first
            % gap has more density than the moment in question then can't
            % interpolate
            
            
            % If the first gap with density and the last are different
            % gaps, interpolate to determine which
            gap_50 = interp1(cumsum(econparams.mu((above_0+1):(reaches_100 + 1))), above_0:reaches_100, 0.5);
            gap_75 = interp1(cumsum(econparams.mu((above_0+1):(reaches_100 + 1))), above_0:reaches_100, 0.75);
            gap_90 = interp1(cumsum(econparams.mu((above_0+1):(reaches_100 + 1))), above_0:reaches_100, 0.9);
            gap_95 = interp1(cumsum(econparams.mu((above_0+1):(reaches_100 + 1))), above_0:reaches_100, 0.95);
            gap_99 = interp1(cumsum(econparams.mu((above_0+1):(reaches_100 + 1))), above_0:reaches_100, 0.99);
        else
            error("Something awry in calculating markup moments")
        end
        
        if isnan(gap_50)
            if econparams.mu(above_0+1) > 0.5
                gap_50 = above_0;
            else
                error("ERROR")
            end
        end
        
        if isnan(gap_75)
            if econparams.mu(above_0+1) > 0.75
                gap_75 = above_0;
            else
                error("ERROR")
            end
        end

        if isnan(gap_90)
            if econparams.mu(above_0+1) > 0.90
                gap_90 = above_0;
            else
                error("ERROR")
            end
        end
        
        if isnan(gap_95)
            if econparams.mu(above_0+1) > 0.95
                gap_95 = above_0;
            else
                error("ERROR")
            end
        end
        
        if isnan(gap_99)
            if econparams.mu(above_0+1) > 0.99
                gap_99 = above_0;
            else
                error("ERROR")
            end
        end

        % If index is Nan then moment is in tied position
        if gap_50 == 0
            econparams.median_gross_markup = econparams.markups(econparams.n + 1);
        elseif gap_50 == econparams.n
            econparams.median_gross_markup = econparams.markups(2*econparams.n + 1);
        % Otherwise, interpolate and weight by distance between markups
        else
            share_near = gap_50 - floor(gap_50);
            econparams.median_gross_markup = (1 - share_near)*econparams.markups(econparams.n + 1 + floor(gap_50)) + ...
                share_near*econparams.markups(econparams.n + 1 + ceil(gap_50));
        end

        % If index is Nan then moment is in tied position
        if gap_90 == 0
            econparams.pct90_gross_markup = econparams.markups(econparams.n + 1);
        elseif gap_90 == econparams.n
            econparams.pct90_gross_markup = econparams.markups(2*econparams.n + 1);
        % Otherwise, interpolate and weight by distance between markups
        else
            share_near = gap_90 - floor(gap_90);
            econparams.pct90_gross_markup = (1 - share_near)*econparams.markups(econparams.n + 1 + floor(gap_90)) + ...
                share_near*econparams.markups(econparams.n + 1 + ceil(gap_90));
        end
        
        % If index is Nan then moment is in tied position
        if gap_75 == 0
            econparams.pct75_gross_markup = econparams.markups(econparams.n + 1);
        elseif gap_75 == econparams.n
            econparams.pct75_gross_markup = econparams.markups(2*econparams.n + 1);
        % Otherwise, interpolate and weight by distance between markups
        else
            share_near = gap_75 - floor(gap_75);
            econparams.pct75_gross_markup = (1 - share_near)*econparams.markups(econparams.n + 1 + floor(gap_75)) + ...
                share_near*econparams.markups(econparams.n + 1 + ceil(gap_75));
        end
        
        % If index is Nan then moment is in tied position
        if gap_95 == 0
            econparams.pct95_gross_markup = econparams.markups(econparams.n + 1);
        elseif gap_95 == econparams.n
            econparams.pct95_gross_markup = econparams.markups(2*econparams.n + 1);
        % Otherwise, interpolate and weight by distance between markups
        else
            share_near = gap_95 - floor(gap_95);
            econparams.pct95_gross_markup = (1 - share_near)*econparams.markups(econparams.n + 1 + floor(gap_95)) + ...
                share_near*econparams.markups(econparams.n + 1 + ceil(gap_95));
        end
        
        % If index is Nan then moment is in tied position
        if gap_99 == 0
            econparams.pct99_gross_markup = econparams.markups(econparams.n + 1);
        elseif gap_99 == econparams.n
            econparams.pct99_gross_markup = econparams.markups(2*econparams.n + 1);
        % Otherwise, interpolate and weight by distance between markups
        else
            share_near = gap_99 - floor(gap_99);
            econparams.pct99_gross_markup = (1 - share_near)*econparams.markups(econparams.n + 1 + floor(gap_99)) + ...
                share_near*econparams.markups(econparams.n + 1 + ceil(gap_99));
        end
    else
        reaches_50 = find(cumsum(econparams.mu) >= 0.5, 1) - 1;
        reaches_75 = find(cumsum(econparams.mu) >= 0.75, 1) - 1;
        reaches_90 = find(cumsum(econparams.mu) >= 0.9, 1) - 1;
        reaches_95 = find(cumsum(econparams.mu) >= 0.95, 1) - 1;
        reaches_99 = find(cumsum(econparams.mu) >= 0.99, 1) - 1;
        econparams.median_gross_markup = econparams.markups(econparams.n + 1 + reaches_50);
        econparams.pct90_gross_markup = econparams.markups(econparams.n + 1 + reaches_90);
    end
        
    % Compute entry rate (AA 2019 Eq 25)
    if econparams.entry_flag ~= -1
        econparams.entry_rate = 0.5*sum(econparams.mu .* flipud(econparams.x_vector_e));
    end
    
    % Compute research demand for each firm
    econparams.research_labor_demand = (econparams.x_vector./econparams.B).^(1/econparams.gamma);
    econparams.research_labor_demand_e = (econparams.x_vector_e./econparams.B_e).^(1/econparams.gamma);
    
    % Compute total labor demand for each firm
    econparams.total_labor_demand = econparams.research_labor_demand + econparams.production_labor_demand;
    
    % Compute total research and production labor demand
    econparams.L_production = sum(econparams.mu.*econparams.production_labor_demand((econparams.n + 1):(2*econparams.n + 1))) + ...
        sum(econparams.mu.*econparams.production_labor_demand((econparams.n + 1):-1:1));
    
    econparams.L_research_noent = sum(econparams.mu.*econparams.research_labor_demand((econparams.n + 1):(2*econparams.n + 1))) + ...
        sum(econparams.mu.*econparams.research_labor_demand((econparams.n + 1):-1:1));
    
    if econparams.entry_flag == -1
        econparams.L_research = econparams.L_research_noent + ...
            econparams.research_labor_demand_e;
    else
        econparams.L_research = econparams.L_research_noent + ...
            sum(econparams.mu.*econparams.research_labor_demand_e((econparams.n + 1):-1:1));
    end
    
    % Compute total labor demand
    econparams.L = econparams.L_production + econparams.L_research;

    % Research productivity
    econparams.research_productivity = econparams.g/econparams.L_research;
    
    % Compute RD to GDP
    econparams.rd_gdp = econparams.omega*econparams.L_research;
    econparams.rd_gdp_noent = econparams.omega*econparams.L_research_noent;
    
    % Compute leader share of RD 
    if econparams.entry_flag == -1
        econparams.rd_leader_share = econparams.research_labor_demand(econparams.n + 1:2*econparams.n + 1)./(econparams.research_labor_demand(econparams.n + 1:2*econparams.n + 1) + ...
            econparams.research_labor_demand(econparams.n + 1:-1:1) + econparams.research_labor_demand_e);
    else
        econparams.rd_leader_share = econparams.research_labor_demand(econparams.n + 1:2*econparams.n + 1)./(econparams.research_labor_demand(econparams.n + 1:2*econparams.n + 1) + ...
            econparams.research_labor_demand(econparams.n + 1:-1:1) + econparams.research_labor_demand_e((econparams.n + 1):-1:1));
    end
    econparams.rd_leader_share(isnan(econparams.rd_leader_share)) = 0;
    econparams.rd_leader_share = sum(econparams.mu .* econparams.rd_leader_share);
    
    % Compute R&D spend to value
    econparams.rd_value = econparams.omega*econparams.research_labor_demand./econparams.value_function;
    
    % Initialize BGP potential output
    econparams.Q0 = 1;
    
    % Initialize BGP consumption
    econparams.Y0 = exp(log(econparams.Q0) - sum(econparams.mu .* log(econparams.markups((econparams.n + 1):(2*econparams.n + 1)))) - log(econparams.omega));
     
    % Compute welfare
    econparams.W = (log(econparams.Y0/econparams.Q0) - econparams.L)/econparams.rho + ...
        (econparams.g/(econparams.rho^2));
    
    % Simulate economy
    if fhk_flag
        econparams = simulate_economy(econparams, false);    
    end
end