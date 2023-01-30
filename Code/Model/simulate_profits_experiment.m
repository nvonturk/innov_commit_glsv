%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: simulate_profits_experiment.m
% Author: Nicholas von Turkovich
% Date: 05/05/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = simulate_profits_experiment(param_file, transition_type)

    %% Gather data

    % Economy variables
    try
        econparams = load("Data/Intermediate/" + param_file + "_econparams.mat");
    catch
        error(strcat("The econparams variable for '", param_file, "' was not found. Ensure that you have fully run the code under switch baseline prior to running this function."))
    end
    
    % Starting (calibrated) BGP
    econparams_bgp = value_function_iteration_robust(econparams, false);

    % Policy to use
    try
        policy = load("Data/Intermediate/" + param_file + "_" + transition_type + ".mat");
    catch
        error(strcat("The policy variable under '", transition_type, "' for the '", param_file, "' case was not found. Ensure that you have fully run the code under switch baseline prior to running this function."))
    end
        
    if transition_type == "consistent"
        policy = policy.heta_optcons_duple;
    elseif transition_type == "commitment"
        policy = policy.heta_opt_duple;
    else
        error("Specify one of two baseline transitions")
    end

    transition = compute_transition(econparams, econparams_bgp, policy, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);

    % Generate draw of events
    seed = 2019;
    rng(seed, 'Threefry')
    N = 500000;
    draws = rand(N, length(transition.heta_path));

    %% Simulate time pattern of proifts
    
    % Innovation comparison occurs at time t and we track profits from t + 1 to
    % t + z
    t_span = [1, 24]./econparams.frequency;
    for t = t_span
        sigma_start = [-5,-4,-3,-2,-1,0,1,2,3,4,5];
        tic
        simulations = simulate_profits(econparams, transition, sigma_start, t + 1, draws);
        toc
        save(strcat("Data/Intermediate/", param_file, "_", transition_type, "_sim_t", num2str(t*econparams.frequency), ".mat"), "-struct", "simulations")
    end
    disp("Completed profits experiment for GLSV")

end



