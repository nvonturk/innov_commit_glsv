%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: generate_econparams.m
% Author: Nicholas von Turkovich
% Date: 11/4/2021
% Note(s): Takes a set of parameters and instantiates the necessary
% data/structures to solve for the BGP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = generate_econparams(params, disp_flag)
    
    % Order of arguments in params:
    %  1.  Frequency (1/(number of periods in year))
    %  2.  Number of gaps (n)
    %  3.  Discount rate (rho, annual decimal)
    %  4.  Productivity multiplier (lambda)
    %  5.  Convexity parameter for R&D costs (gamma)
    %  6.  Elasticity of substitution across varieties (kappa)
    %  7.  Rate of patent expiry (heta, annual decimal)
    %  8.  Fraction of time patent expiry wipes out lead vs. one step (zeta, decimal share)
    %  9.  Tax rate on profits (tau, decimal share)
    %  10. Tax subsidy on R&D (tau_rd, decimal share)
    %  11. R&D cost scaling parameter for incumbents (B, annual)
    %  12. R&D cost scaling parameter for entrants (B_e, annual)
    %  13. Rate of quick catch up after innovation (phi, decimal share)
    %  14. Rate of quick catch up for entrants after innovation (phi_e, decimal share)
    %  15. Leapfrogging for incumbents (leap)
    %  16. Leapfrogging for entrants (leap_e)
    %  17. Maximum (initial) innovation intensity (x_bar)
    %  18. Entry flag (1 = entry, 0 = no entry)
    %  19. Elastic labor flag (1 = elastic, 0 = inelastic (omega \neq 1))
    %  20. Tolerance on value function difference for VFI
    %  21. Tolerance on wage rate for VFI with inelastic labor supply
    
    %% Setup econparams object
    
    % From the inputs
    econparams.frequency = params(1);
    econparams.n = params(2);
    econparams.rho = power(1 + params(3), econparams.frequency) - 1;
    econparams.lambda = params(4);
    econparams.gamma = params(5);
    econparams.kappa = params(6);
    econparams.heta = params(7)*econparams.frequency;
    econparams.zeta = params(8);
    econparams.tau = params(9);
    econparams.tau_rd = params(10);
    econparams.B = params(11)*econparams.frequency;
    econparams.B_e = params(12)*econparams.frequency;
    econparams.phi = params(13);
    econparams.phi_e = params(14);
    econparams.leap = params(15);
    econparams.leap_e = params(16);
    econparams.x_bar = params(17);
    econparams.entry_flag = params(18);
    econparams.elastic_labor_flag = params(19);
    econparams.v_tol = params(20);
    econparams.w_tol = params(21);
    econparams.disp_flag = disp_flag;
    econparams.params = params;
    
    %% Initial values for objects
    
    % Value function where index = 1 is laggard in largest gap, n+1 is
    % tied, and 2n + 1 is leader in largest gap
    econparams.value_function = ones((econparams.n*2 + 1),1);
    
    % Same structure as value function but for innovation intensities
    econparams.x_vector = zeros((econparams.n*2 + 1), 1);
    
    if econparams.entry_flag == -1
        % Undirected entry
        econparams.x_vector_e = 0;
    else
        % Potential entrant for each industry gap
        econparams.x_vector_e = zeros((econparams.n + 1), 1);
    end
    
    % If undirected entry, need to instantiate a mu
    if econparams.entry_flag == -1
        econparams.mu = (1/(econparams.n+1))*ones(econparams.n + 1, 1);
    end

    % Wage is 1 in equilibrium for case of inelastic labor supply
    if ~econparams.elastic_labor_flag
        econparams.w = 1;
        econparams.omega = 0.9;
    else
        econparams.omega = 1;
    end
    
    %% Set up transition matrices
    econparams = generate_transition_matrices(econparams);
    
    %% Calculate profits
    
    % If kappa is -1, flag to have fully imperfect substitutes
    if econparams.kappa ~= -1
    
        % Vector of relative prices for each industry gap
        relative_prices = -1*ones((econparams.n+1),1);
        for i = 1:(econparams.n+1)
            productivity_gap = i - 1;
            relative_prices(i) = fzero(@(v_guess) ((1/econparams.lambda)^productivity_gap)*...
            (econparams.kappa*v_guess^(econparams.kappa - 1) + 1)/...
            (econparams.kappa + v_guess^(econparams.kappa - 1)) - ...
            v_guess^econparams.kappa, 1);
        end

        % Vector of markups
        econparams.markups = zeros((econparams.n*2 + 1),1);
        for i = 0:econparams.n
           econparams.markups(econparams.n + 1 + i) = (econparams.kappa + relative_prices(i+1)^(1 - econparams.kappa))/...
               (econparams.kappa - 1);
           econparams.markups(econparams.n + 1 - i) = (econparams.kappa*relative_prices(i+1)^(1 - econparams.kappa) + 1)/...
               ((econparams.kappa - 1)*(relative_prices(i+1)^(1 - econparams.kappa)));
        end

        % Vector of profits
        econparams.profits = 1 - econparams.kappa/(econparams.kappa - 1)*(1./econparams.markups);
        
    else
        
        % Vectors of markups and profits
        econparams.markups = ones((econparams.n*2 + 1),1);
        econparams.profits = zeros((econparams.n*2 + 1),1);
        for s = 0:econparams.n
            econparams.markups(econparams.n + 1 + s) = econparams.lambda^s;
            econparams.profits(econparams.n + 1 + s) = 1 - econparams.lambda^(-s);
        end
        
    end
    
    % Vector of labor demand
    econparams = compute_labor_demand(econparams);

end