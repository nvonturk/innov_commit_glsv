%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: plot_core.m
% Author: Nicholas von Turkovich
% Date: 7/6/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_core(param_file)
    
    %% Script setup

    % Set up specification for charts
    specification = param_file;

    % Set up where to pull data from
    directory_data = "Data/Intermediate/";
    filestart_data = directory_data + specification;

    % Set up where to place figures
    directory_fig = "Output/Figures/Main/";
    filestart_fig = directory_fig + specification;

    %% Gather data

    % Load values from running policy exercises

    % Economy variables
    econparams = load(filestart_data + "_econparams.mat"); % Econparams
    econparams_bgp = load(filestart_data + "_econparams_bgp.mat"); % Calibrated BGP

    % Policies with a single reoptimization
    commitment = load(filestart_data + "_commitment.mat"); % Commitment policies
    contour_mat = load(filestart_data + "_contour.mat"); % Contour plot with commitment duples
    consistent = load(filestart_data + "_consistent.mat"); % Consistent policies

    fprintf("Best response to optimal commitment policy of second social planner = %1.2f%% \n", commitment.sp2_best_response*100/econparams.frequency)

    % Policies with two reoptimizations
    heta_triple_commitment = load(filestart_data + "_heta_triple_commitment.mat").heta_triple_commitment; % Optimal triple with commitment
    tripcons = load(filestart_data + "_tripcons.mat"); % Consistent policies (double reoptimization)

    % Transitions and subsequent data with a single reoptimization

    % Given policies loaded, compute the transitions for plotting
    duple_commitment_transition = compute_transition(econparams, econparams_bgp, commitment.heta_duple_commitment, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);
    optcons_duple_transition = compute_transition(econparams, econparams_bgp, consistent.heta_optcons_duple, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);

    % CEQ gain from moving from optimal consistent to optimal commitment
    % (data point in text)
    optcons_commit_ceq = compute_ceq_change(optcons_duple_transition.W, duple_commitment_transition.W, econparams.rho);
    fprintf("CEQ gain with single reoptimization moving from consistent to commitment policy = %1.2f %%\n", optcons_commit_ceq*100)

    % Given policies loaded, compute terminal BGPs
    duple_commitment_terminal = value_function_iteration_robust(update_heta(econparams, commitment.heta_duple_commitment(2)), false);
    optcons_duple_terminal = value_function_iteration_robust(update_heta(econparams, consistent.heta_optcons_duple(2)), false);

    % Transitions and subsequent data with two reoptimizations

    % Given policies loaded, compute the transitions for plotting
    triple_commitment_transition = compute_transition(econparams, econparams_bgp, heta_triple_commitment, repmat(econparams.rho,1,3), econparams.interval_lengths(end-2:end));
    optcons_triple_transition = compute_transition(econparams, econparams_bgp, tripcons.heta_max_optcons, repmat(econparams.rho,1,3), econparams.interval_lengths(end-2:end));

    % CEQ gain from moving from optimal consistent to optimal commitment
    % (data point in paper)
    optcons_commit_triple_ceq = compute_ceq_change(optcons_triple_transition.W, triple_commitment_transition.W, econparams.rho);
    fprintf("CEQ gain with double reoptimization moving from consistent to commitment policy = %1.2f %%\n", optcons_commit_triple_ceq*100)

    % Given policies loaded, compute terminal BGPs
    triple_commitment_terminal = value_function_iteration_robust(update_heta(econparams, heta_triple_commitment(3)), false);
    optcons_triple_terminal = value_function_iteration_robust(update_heta(econparams, tripcons.heta_max_optcons(3)), false);

    % Simulations for the time pattern of profits exercise
    sims = {load(filestart_data + "_consistent_sim_t1.mat"),...
        load(filestart_data + "_consistent_sim_t24.mat")};

    % Output from the strategic decomp exercise
    strat_decomp = load(filestart_data + "_consistent_stratdecomp.mat");

    % State dependent policy exercises
    commitment_statedependent = load(filestart_data + "_hub40_commitment_statedep1.mat");
    consistent_statedependent = load(filestart_data + "_hub40_consistent_statedep1.mat");

    % Data structures for the pre-sclerosis calibration
    specification_presclerosis = specification + "_presclerosis";
    filestart_data_presclerosis = directory_data + specification_presclerosis;
    econparams_presclerosis = load(filestart_data_presclerosis + "_econparams.mat"); 
    commitment_presclerosis = load(filestart_data_presclerosis + "_commitment.mat");
    contour_mat_presclerosis = load(filestart_data_presclerosis + "_contour.mat"); 
    consistent_presclerosis = load(filestart_data_presclerosis + "_consistent.mat"); 

    %% Tables

    % Calibration table and policy path tables
    tables(specification, econparams, optcons_duple_transition, optcons_duple_terminal, duple_commitment_transition, duple_commitment_terminal,...
        optcons_triple_transition, optcons_triple_terminal, triple_commitment_transition, triple_commitment_terminal, econparams_presclerosis, ...
        commitment_statedependent, consistent_statedependent)

    %% Figures

    % Figures organized by appearance in paper

    % "Main" plot which details the the difference in commitment and consistent
    % policy in (eta1, eta2) space with the added dimension of welfare in a
    % contour plot
    figure_isowelfare_contours(econparams, commitment, contour_mat, consistent, specification);

    % Visual illustration of why there exists a (unique) consistent policy for
    % each choice of eta1
    figure_sp2_best_response(econparams, contour_mat, consistent)

    % Illustration of mu_T for different time consistent duples
    figure_mu_transition(econparams, econparams_bgp, consistent)

    % Series of plots of main equilibrium objects at different points in the
    % commitment/consistent transitions
    figure_menu_xmu(econparams, duple_commitment_transition, ...
        optcons_duple_transition)

    % Time series of growth/markups for the optimal time-consistent transition
    figure_optimal_consistent(econparams, econparams_bgp, optcons_duple_transition);

    % Research productivity
    figure_research_productivity(econparams, econparams_bgp, duple_commitment_transition, optcons_duple_transition);

    % Changes in strategies given increase in future (long-run) patent expiry
    % rate; changes in welfare associated with an increase in strategies
    figure_time_consistency_components(econparams, econparams_bgp, consistent, optcons_duple_transition);

    % Strategic decomposition of transition
    figure_strategic_decomp(strat_decomp);

    % Time pattern of profits following an innovation at a particular
    % technology position
    figure_time_pattern_profits(sims);

    % Core isowelfare plot adjusting a single parameter at a time
    figure_isowelfare_contours_robust(filestart_data, econparams, commitment, consistent);

    % Plots detailing the effect of shifting T
    figure_adjusting_t(filestart_data, econparams)

    % Produces a 2x1 plot with consistent/commitment policies conditional on a
    % choice of eta1
    figure_triple_consistent(econparams, tripcons, heta_triple_commitment)

    % Figures in the IA

    % Illustrates the accuracy of using a spline to fill in consistent duple
    % values
    figure_spline_approx(econparams, consistent);

    % Illustrates loss associated with various triples and the response of the
    % second social planner
    figure_triple_response(specification)

    % "Main" plot but for the pre-sclerosis economy
    figure_isowelfare_contours(econparams_presclerosis, commitment_presclerosis, contour_mat_presclerosis, consistent_presclerosis, specification_presclerosis);

    % Plot shows decomposition of welfare based on markup
    % distortions/disutility of labor and productivity growth
    figure_consistent_tangent(econparams, commitment, consistent);

    % State dependency under consistent policy
    figure_state_dependency(econparams, commitment_statedependent, consistent_statedependent);

    % Core plot adjusting T
    figure_isowelfare_contours_adjusting_t(filestart_data, econparams, commitment, consistent)

end
